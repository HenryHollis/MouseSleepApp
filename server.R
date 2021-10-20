#check about DeltaF output format

#assumptions: -last 3 rows of all xlsx files are throwaway
#             -wake files have w or W in them


library(shiny)

# use the below options code if you wish to increase the file input limit
options(shiny.maxRequestSize = 100*1024^2)

shinyServer(function(input,output, session) {
  
  
  gather_accross_files = function(){
    # we want to get cell statistics across ALL the files uploaded...
    wake_file_idx = grep("[wW]", v, perl = T, value = F)
    #validate(need(length(wake_file_idx)!=0, "Please upload wake files (has w or W in filename)" )) 
    len = length(input$file$name)
    data =  read_excel(path = input$file$datapath[1], skip = 1, trim_ws = TRUE)
    ncols = length(data)
    updateSliderInput(session, "num", max = ncols-1)
    n = dim(data)[1]
    data = data[1:(n-3),-1]  #last three rows of data are always throwaways
    data = drop_na(data)
    max_amp = zeros(1,ncols)
    if (length(wake_file_idx)==0 | 1 %in% wake_file_idx){  #if no wake files get max amp from first file
      local_medians = apply(data, 2, median)  #apply median to every cell col in data
      local_maxs = apply(data, 2, max)
      local_amp = local_maxs-local_medians
      replace = local_amp > max_amp
      max_amp[replace] = local_amp[replace]
    }
    
    if (len > 1){
      for (i in 2:len){
        local_data = read_excel(path = input$file$datapath[i], skip = 1, trim_ws = TRUE)
        validate(need(ncols == length(local_data), "Please Select Files with the same number of columns" ))  #make sure data has equal cols
        ncols = length(local_data)
        local_data = local_data[,-1]
        local_data = drop_na(local_data)
        
        if (i %in% wake_file_idx){
          local_medians = apply(data, 2, median)  #apply median to every cell col in data
          local_maxs = apply(data, 2, max)
          local_amp = local_maxs-local_medians
          replace = local_amp > max_amp
          max_amp[replace] = local_amp[replace]
        }
        
        rbind(data, local_data)
      }
    }
    
    cell_mins = apply(data, 2, min)

    data = mapply("-", data, cell_mins )
    data = as.data.frame(data)
    cell_maxs = apply(data, 2, max)
    data = mapply("/", data, cell_maxs)
    
    cell_baseline = apply(data, 2, median)
    cell_variance = apply(data, 2, mad)

    return(rbind(cell_mins, cell_maxs, cell_baseline, cell_variance, max_amp))
  }
  
  
  get_results = function(data) {

    # method for doing the analysis on an xlsx file
    threshold = input$thresholdFactor
    data = drop_na(data)
    Time = unname(unlist(data[,1]))
    
    cells = select(data, -1)
    cell_stats = gather_accross_files()  #get min, max, baseline and var from all files
    get_area = function(col, stats){
      #my_min = stats[1]
      #col_norm = col - my_min
      max_amp = stats[5]
      col = col / max_amp
      baseline = stats[3]   # get median of all cell data
      variance = stats[4]    
      
      all_spike_locations = (col > input$thresholdFactor*variance + baseline)

      ## method from https://masterr.org/r/how-to-find-consecutive-repeats-in-r/
      rle.all_spike_locations = rle(as.numeric(all_spike_locations))
      myruns = which(rle.all_spike_locations$values == TRUE & rle.all_spike_locations$lengths >= input$lengthCutoff)
      runs.lengths.cumsum = cumsum(rle.all_spike_locations$lengths)
      ends = runs.lengths.cumsum[myruns]
      newindex = ifelse(myruns>1, myruns-1, 0)
      starts = runs.lengths.cumsum[newindex] + 1   #starts indices of runs > 10 in length
      
      if (0 %in% newindex) starts = c(1,starts)
      
      
      cum_area = 0
      if (length(ends)==0) {
        return(c(0.0,0.0,0.0,0.0,0.0, 0.0,input$thresholdFactor*variance+baseline, 0.0))
        
      }else{
        max_spike_area = 0.0
        longest_time = 0.0
        #lin_coeff = c()
        exp_coeff = c()
        all_falling_spike_duration = c()
        all_rising_spike_duration = c()
        
        for(i in 1:length(ends)){                      # for every spike...
          time = Time[starts[i]:ends[i]]               # time interval for spike_i
          
          if (Time[ends[i]]-Time[starts[i]] > longest_time){   # identify longest spike-time interval
            longest_time = Time[ends[i]]-Time[starts[i]]
          }
          
          spike_points = col[starts[i]:ends[i]] # y values of spikes
          
          max_spike_val = max(spike_points)
          max_spike_idx = which(spike_points == max_spike_val) + starts[i]
          
          falling_spike_times = Time[max_spike_idx:ends[i]]               # spike time-interval, times when spike is declining
          falling_spike_vals = col[max_spike_idx:ends[i]]  # values of spike as it's declining
          rising_spike_times = Time[starts[i]:max_spike_idx]
          
          rising_spike_duration = range(rising_spike_times)[2] - range(rising_spike_times)[1]
          falling_spike_duration = range(falling_spike_times)[2] - range(falling_spike_times)[1]
          all_falling_spike_duration = cbind(all_falling_spike_duration, falling_spike_duration)
          all_rising_spike_duration = cbind(all_rising_spike_duration, rising_spike_duration)
        
         # lm_result = lm(falling_spike_vals ~ falling_spike_times)        #run linear regression

          exp_result = lm(log(falling_spike_vals) ~ falling_spike_times)  #run "exp regression"
          
          #lin_coeff = cbind(lin_coeff, unname(lm_result$coefficients[2]))
          exp_coeff = cbind(exp_coeff, unname(exp_result$coefficients[2]))
          
          
          spike_area = trapz(time, spike_points)    # integrate spike with trapz package (trapazoidal integration)
          
          if (spike_area > max_spike_area){          # identify the max spike area
            max_spike_area = spike_area
          }
          cum_area = cum_area+spike_area
        }
        #mean_lin_coeff = mean(lin_coeff, na.rm=TRUE)
        mean_exp_coeff = mean(exp_coeff, na.rm=TRUE)
        mean_rising_spike_duration = mean(all_rising_spike_duration, na.rm = T)
        mean_falling_spike_duration = mean(all_falling_spike_duration, na.rm = T)
        
        
        results = c(  
          cum_area/length(ends),                        # average spike area
          length(ends),                                 # number of spikes
          max_spike_area,                               # max_spike_area
          longest_time,                                 # longest time
          mean_rising_spike_duration,
          mean_falling_spike_duration,
          input$thresholdFactor*variance+baseline,      # threshold
          mean_exp_coeff                               # average coeff from exp linear regression
         # mean_lin_coeff                                # avg coeff from simple linear regression
          
        )
        return(results)
      }
    }
    
    results = mapply(get_area, as.data.frame(cells), as.data.frame(cell_stats))
    results = t(results)
    results = as_tibble(results) %>% mutate(cell_id = row_number()) %>% select( cell_id, everything())
    colnames(results) = c("cell_id", "avg_spike_area", "num_spikes", "max_spike_area", "longest_spike", "mean_rising_spike_duration", "mean_falling_spike_duration","threshold", "mean_exp_coeff")
    
    return(as_tibble(results))
    
  }

  
  get_Fdelta = function(data){
    cell_stats = gather_accross_files()
    data = drop_na(data)
    cells = select(data, -1)
    Time = unname(unlist(data[,1]))
    
    
    helper = function(col, stats){
      #my_min = stats[1]
      #col = col - my_min
      max_amp = stats[5]
      col = col / max_amp
      
      median = median(col)
      col = (col - median)/median
      
      return(col)
    }
    cell_out = mapply(helper, as.data.frame(cells), as.data.frame(cell_stats))
    #cell_out = as_tibble(cell_out) %>% mutate(Time = data[,1]) %>% select( Time, everything())
    #cell_out = t(cell_out)
  }
  
  get_spike_inx = function(col, baseline, var){
    #for creating a logical vector of where spikes are so that they can be displayed on plot
    all_spike_locations = (col > input$thresholdFactor*var + baseline)

    ## method from https://masterr.org/r/how-to-find-consecutive-repeats-in-r/
    rle.all_spike_locations = rle(as.numeric(all_spike_locations))          
    for (i in 1:length(rle.all_spike_locations$lengths)){
      if(rle.all_spike_locations$lengths[i] < input$lengthCutoff){          #TODO FIXME
        rle.all_spike_locations$values[i] = 0
      }
    }
    
    return(inverse.rle(rle.all_spike_locations))
    
  }
  
  ## input$file is a data frame and contains the details around the name, size and temp location of the files uploaded
  # this reactive output display the content of the input$file dataframe
  output$filedf <- renderTable({
    if(is.null(input$file)){return ()}
    input$file # the file input data frame object that contains the file attributes
  })
  
  # Extract the file path for file
  output$filedf2 <- renderTable({
    if(is.null(input$file)){return ()}
    input$file$datapath # the file input data frame object that contains the file attributes
  })
  
  ## Below code to display the structure of the input file object
  output$fileob <- renderPrint({
    if(is.null(input$file)){return ()}
    str(input$file)
  })
  
  ## Side bar select input widget coming through renderUI()
  # Following code displays the select input widget with the list of file loaded by the user
  output$selectfile <- renderUI({
    if(is.null(input$file)) {return()}
    
    list(hr(), 
         helpText("Select the files for which you need to see data and summary stats"),
         selectInput("Select", "Select", choices=input$file$name)
    )
    
  })
  
  
  
  
  
  
  #run the get_results function on all the files uploaded for downloading later....
  run_all = function(){
    files = NULL;
    withProgress(message = 'Processing data...', value = 0, {
      
    #loop through the sheets
    len = length(input$file$name)

    for (i in 1:len){
      #write each sheet to a csv file, save the name
      #trunc_name = substr(input$file$name[i], 1, nchar(input$file$name)-5)
      trunc_name = rmv.ext(input$file$name[i])
      fileName1 = paste(trunc_name,"spike_stats.csv",sep = "_")
      fileName2 = paste(trunc_name,"deltaF.csv",sep = "_")

      data = read_excel(path = input$file$datapath[i], skip = 1, trim_ws = TRUE)
      temp = get_results(data)
      temp2 = (get_Fdelta(data))
      write.table(temp, fileName1, sep = ',', row.names = F, col.names = T)
      write.table(temp2, fileName2, sep = ',', row.names = T, col.names = F)
      
      files = c(fileName1,files)
      files = c(fileName2, files)
      incProgress(1/len, detail = paste("Processing file", i))
      
    }
  })
    return(files)
  }
  
  ## Summary Stats code ##
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$results <- renderPrint({
    if(is.null(input$file)){return()}
    as.data.frame(get_results(read_excel(path = input$file$datapath[input$file$name==input$Select], skip = 1, trim_ws = TRUE)))
    

  })

  
  get_data <- reactive({
    if(is.null(input$file)){return()}
    data = read_excel(path = input$file$datapath[input$file$name==input$Select], skip = 1, trim_ws = TRUE)
    data = drop_na(data)
    #data[,-1] = data[,-1]-min(data[,-1])
    
    return(data)
  })
  
  
  # function for plotting
  output$plot <- renderPlot({ 
    
    cell_stats = gather_accross_files()  #get min, max, baseline and var from all files
    data = get_data()
    Time = data[, 1]
    data = data[, -1]
    
    cell = unname(unlist(data[ ,input$num]))   #columns of data, user selected
    #my_min = cell_stats[1, input$num]
    #cell_norm = cell - my_min
    max_amp = cell_stats[5, input$num]
    cell = cell / max_amp
    
    baseline = cell_stats[3, input$num]
    variance = cell_stats[4, input$num]
    
    logical = as.logical(get_spike_inx(cell, baseline, variance))

    
    ggplot(data, aes(x = as.matrix(Time), y = cell))+
      geom_point(aes(color = logical))+
      scale_color_manual(labels = c("<baseline", ">baseline"), values=c('blue', 'red')) +
      labs(color = "Normalized Spikes")+
      xlab("Time")
    
    # ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) + 
    #     geom_point(mapping = aes(color = drv)) + 
    
    
  })
  
  ## MainPanel tabset renderUI code ##
  # the following renderUI is used to dynamically generate the tabsets when the file is loaded. 
  # Until the file is loaded, app will not show the tabset.
  output$tb <- renderUI({
    if(is.null(input$file)) {return()}
    else
      tabsetPanel(
        tabPanel("Input File Object DF ", tableOutput("filedf"), tableOutput("filedf2")),
        tabPanel("Plot", plotOutput("plot")),
        tabPanel("File Results", verbatimTextOutput("results")),
        tabPanel("Input File Object Structure", verbatimTextOutput("fileob"))
        
      )
  })
  
  #Code to zip together the analyses into a .zip archive and download them
  output$download <- downloadHandler(
    filename = function(){
      paste0("analysis.zip")
      
    },
    content = function(file){
      #go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files = run_all()
      
      #create the zip file
      withProgress(message = 'Zipping files', value = .5, {
      zip(file,files)
      })
    }
  )
  
  
})