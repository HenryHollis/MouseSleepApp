
library(shiny)

# use the below options code if you wish to increase the file input limit
options(shiny.maxRequestSize = 100*1024^2)

shinyServer(function(input,output, session) {

  
  rm_cells = function(string){
    rmcells = str_replace_all(string, " ", "")  #cleans string specifying which cells to remove
    rmcells = as.numeric(strsplit(rmcells, ",")[[1]])
    rmcells = unique(rmcells)
  }
  

  gather_across_files = function(){
    max_file_idx = grep("[Mm][Aa][xX]", input$file$name, perl = T, value = F)
    validate(need(length(max_file_idx)==1, "Please upload a single file with max values for each cell" ))  #make sure data has equal cols
    files = input$file$name[-max_file_idx]
    usable_input_file_idx = which(input$file$name %in% files)
    len = length(usable_input_file_idx)
    data =  read_excel(path = input$file$datapath[ usable_input_file_idx[1]], skip = 1, trim_ws = TRUE)  # read the first file
    ncols = length(data)
    updateSliderInput(session, "num", max = ncols-1)
    
    if (len > 1){        #if there are multiple files being uploaded:
      for (i in usable_input_file_idx[-1]){
          local_data = read_excel(path = input$file$datapath[i], skip = 1, trim_ws = TRUE)
          validate(need(ncols == length(local_data), "Please Select Files with the same number of columns" ))  #make sure data has equal cols
          ncols = length(local_data)
        
      }
    }
    max = read_excel(input$file$datapath[max_file_idx])
    cell_maxs = as.numeric(max)
    validate(need(ncols-1 == length(cell_maxs), "Please make sure max file has same #cells as all data files" ))  #make sure data has equal cols
    validate(need(is.numeric(cell_maxs), "Please make sure max file is numeric" ))
    return(cell_maxs)
  }
  

  get_results = function(data) {
    # method for doing the analysis on an xlsx file
    threshold = input$thresholdFactor
    first_NA_row = min(which(is.na(data), arr.ind = T)[, 1])  # first row with NA values, want to remove all rows after this one
    if(is.finite(first_NA_row)){
      data = data[1:(first_NA_row-1),]                        # removes all info after the first NA, needed for the sample files I got to test this on
    }
    data = drop_na(data)
    Time = unname(unlist(data[,1]))
    
    cells = select(data, -1)
    rmcells = rm_cells(input$rmcells)
    cells[,rmcells] = 0.0
    
    maxs = gather_across_files()  #get min, max, baseline and var from all files
    
    get_area = function(col, maxs){
      baseline = median(col)
      variance = mad(col)
      threshold_line = threshold*variance + baseline

      col = col - threshold_line
      my_max = maxs   #comes from max file data
      col = col / my_max
      
      
      col_for_exp = col - median(col)
      all_spike_locations = (col > 0)

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
        return(c(0.0,0.0,0.0,0.0,0.0,0.0, 0.0,threshold_line, 0.0))
        
      }else{
        max_spike_area = 0.0
        longest_time = 0.0
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
          max_spike_idx = which(spike_points == max_spike_val) + starts[i]  #argmax of max_spike_val
          
          falling_spike_times = Time[max_spike_idx:ends[i]]               # spike time-interval, times when spike is declining
          falling_spike_vals = col_for_exp[max_spike_idx:ends[i]]                 # values of spike as it's declining
          rising_spike_times = Time[starts[i]:max_spike_idx]
          
          rising_spike_duration = range(rising_spike_times)[2] - range(rising_spike_times)[1]
          falling_spike_duration = range(falling_spike_times)[2] - range(falling_spike_times)[1]
          all_falling_spike_duration = cbind(all_falling_spike_duration, falling_spike_duration)
          all_rising_spike_duration = cbind(all_rising_spike_duration, rising_spike_duration)
          
          
          exp_decay_bool = sum(diff(falling_spike_vals)) <2 #if falling spike has 2 or more increases, dont calc decay
          if(exp_decay_bool){
              exp_result = lm(log(falling_spike_vals) ~ falling_spike_times)  #run "exp regression"
              exp_coeff = cbind(exp_coeff, unname(exp_result$coefficients[2]))
          }else{
              exp_coeff = cbind(exp_coeff, NA)
          }
          
          spike_area = trapz(time, spike_points)    # integrate spike with trapz package (trapazoidal integration)
          
          if (spike_area > max_spike_area){          # identify the max spike area
            max_spike_area = spike_area
          }
          cum_area = cum_area+spike_area
        }
        mean_exp_coeff = mean(exp_coeff, na.rm=TRUE)
        mean_rising_spike_duration = mean(all_rising_spike_duration, na.rm = T)
        mean_falling_spike_duration = mean(all_falling_spike_duration, na.rm = T)
        
        
        results = c(  
          cum_area/length(ends),                        # average spike area
          cum_area/(length(ends) * 0.05 * length(col)), #average area / 0.05 * number of rows (accounts for different recording times)
          length(ends),                                 # number of spikes
          max_spike_area,                               # max_spike_area
          longest_time,                                 # longest time
          mean_rising_spike_duration,
          mean_falling_spike_duration,
          input$thresholdFactor*variance+baseline,      # threshold
          mean_exp_coeff                               # average coeff from exp linear regression
        )
        return(results)
      }
    }
    
    results = mapply(get_area, as.data.frame(cells), maxs)
    results = t(results)
    results = as_tibble(results) %>% mutate(cell_id = row_number()) %>% select( cell_id, everything())
    colnames(results) = c("cell_id", "avg_spike_area", "corrected_apike_area", "num_spikes", "max_spike_area", "longest_spike", "mean_rising_spike_duration", "mean_falling_spike_duration","threshold", "mean_exp_coeff")
    
    return(as_tibble(results))
    
  }

  
  get_Fdelta = function(data){
    maxs = gather_across_files()
    first_NA_row = min(which(is.na(data), arr.ind = T)[, 1])  # first row with NA values, want to remove all rows after this one
    if(is.finite(first_NA_row)){
      data = data[1:(first_NA_row-1),]                        # removes all info after the first NA, needed for the sample files I got to test this on
    }
    data = drop_na(data)
    cells = select(data, -1)
    rmcells = rm_cells(input$rmcells)
    cells[,rmcells] =0.0      #set cells labeled as rmcells to 0
    Time = unname(unlist(data[,1]))
    
    
    helper = function(col, maxs){
      baseline = median(col)
      variance = mad(col)
      threshold_line = input$thresholdFactor*variance + baseline
      col = col - threshold_line
      my_max = maxs   #comes from max file data
      col = col / my_max
      
      return(col)
    }
    cell_out = mapply(helper, as.data.frame(cells), maxs)
    colnames(cell_out) = paste( "Cell", seq(1,dim(cell_out)[2]))
    return(cell_out)
  }
  
  get_spike_inx = function(col, baseline, var){
    #for creating a logical vector of where spikes are so that they can be displayed on plot
    all_spike_locations = (col > 0)

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
    max_file_idx = grep("[Mm][Aa][xX]", input$file$name, perl = T, value = F)
    for (i in 1:len){
      if(input$file$datapath[i] != input$file$datapath[ max_file_idx] ){
      #write each sheet to a csv file, save the name
      trunc_name = rmv.ext(input$file$name[i])
      fileName1 = paste(trunc_name,"spike_stats.csv",sep = "_")
      fileName2 = paste(trunc_name,"deltaF.csv",sep = "_")

      data = read_excel(path = input$file$datapath[i], skip = 1, trim_ws = TRUE)
      temp = get_results(data)
      temp2 = (get_Fdelta(data))
      write.table(temp, fileName1, sep = ',', row.names = F, col.names = T)
      write.table(temp2, fileName2, sep = ',', row.names = T, col.names = NA)
      
      files = c(fileName1,files)
      files = c(fileName2, files)
      incProgress(1/len, detail = paste("Processing file", i))
      }
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
    first_NA_row = min(which(is.na(data), arr.ind = T)[, 1])  # first row with NA values, want to remove all rows after this one
    if(is.finite(first_NA_row)){
      data = data[1:(first_NA_row-1),]                        # removes all info after the first NA, needed for the sample files I got to test this on
    }
    data = drop_na(data)

    return(data)
  })
  
  
  # function for plotting
  output$plot <- renderPlot({ 
    
    maxs = gather_across_files()  #get min, max, baseline and var from all files
    data = get_data()
    Time = data[, 1]
    data = data[, -1]

    cell = unname(unlist(data[ ,input$num]))   #columns of data, user selected
    baseline = median(cell)
    variance = mad(cell)
    
    threshold_line = input$thresholdFactor*variance + baseline
    cell = cell - threshold_line
    my_max = maxs[input$num]
    cell = cell / my_max

    logical = as.logical(get_spike_inx(cell, baseline, variance))

    
    ggplot(data, aes(x = as.matrix(Time), y = cell))+
      geom_point(aes(color = logical))+
      scale_color_manual(labels = c("<baseline", ">baseline"), values=c('blue', 'red')) +
      labs(color = "Normalized Spikes")+
      xlab("Time")

    
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