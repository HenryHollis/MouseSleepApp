library(shiny)

# use the below options code if you wish to increase the file input limit
options(shiny.maxRequestSize = 30*1024^2)

shinyServer(function(input,output) {
  
  get_results = function(data) {
    # method for doing the analysis on an xlsx file
    threshold = input$thresholdFactor
    data = drop_na(data)
    #data[,-1] = data[,-1]-min(data[,-1])
    Time = unname(unlist(data[,1]))
    
    cells = select(data, -1)
    
    
    get_area = function(col){
      my_min = min(col)
      col_norm = col - my_min
      my_max = max(col_norm)
      col = col_norm / my_max
      baseline = median(col)   # get median of all cell data
      variance = mad(col)      # measure of variance (adds 1.4826 factor)
      # col_minus_baseline = (col - baseline)
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
        return(c(0.0,0.0,0.0,0.0,input$thresholdFactor*variance+baseline, 0.0, 0.0))
        
      }else{
        max_spike_area = 0.0
        longest_time = 0.0
        lin_coeff = c()
        exp_coeff = c()
        
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
          
          lm_result = lm(falling_spike_vals ~ falling_spike_times)        #run linear regression
          
          exp_result = lm(log(falling_spike_vals) ~ falling_spike_times)  #run "exp regression"
          
          lin_coeff = cbind(lin_coeff, unname(lm_result$coefficients[2]))
          exp_coeff = cbind(exp_coeff, unname(exp_result$coefficients[2]))
          
          
          spike_area = trapz(time, spike_points)    # integrate spike with trapz package (trapazoidal integration)
          
          if (spike_area > max_spike_area){          # identify the max spike area
            max_spike_area = spike_area
          }
          cum_area = cum_area+spike_area
        }
        mean_lin_coeff = mean(lin_coeff, na.rm=TRUE)
        mean_exp_coeff = mean(exp_coeff, na.rm=TRUE)
        
        
        results = c(  
          cum_area/length(ends),                        # average spike area
          length(ends),                                 # number of spikes
          max_spike_area,                               # max_spike_area
          longest_time,                                 # longest time
          input$thresholdFactor*variance+baseline,      # threshold
          mean_exp_coeff,                               # average coeff from exp linear regression
          mean_lin_coeff                                # avg coeff from simple linear regression
          
        )
        return(results)
      }
    }
    
    results = apply(cells, 2, get_area)
    results = t(results)
    results = as_tibble(results) %>% mutate(cell_id = row_number()) %>% select( cell_id, everything())
    colnames(results) = c("cell_id", "avg_spike_area", "num_spikes", "max_spike", "longest_spike", "threshold", "mean_lin_coeff", "mean_exp_coeff")
    
    return(as_tibble(results))
    
  }
  
  get_spike_inx = function(col){
    #for creating a logical vector of where spikes are so that they can be displayed on plot
    
    baseline = median(unlist(col))   # get median of all cell data
    variance = mad(unlist(col))      # measure of variance (adds 1.4826 factor)
    #col_minus_baseline = (col - baseline)
    all_spike_locations = (col > input$thresholdFactor*variance + baseline)
    
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
    
    #loop through the sheets
    for (i in 1:length(input$file$name)){
      #write each sheet to a csv file, save the name
      trunc_name = substr(input$file$name, 1, nchar(input$file$name)-5)
      fileName = paste(trunc_name,".csv",sep = "")
      temp = get_results(read_excel(path = input$file$datapath[i], skip = 1, trim_ws = TRUE))
      write.table(temp, fileName, sep = ',', row.names = F, col.names = T)
      
      files = c(fileName,files)
    }
    return(files)
  }
  
  ## Summary Stats code ##
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$results <- renderPrint({
    if(is.null(input$file)){return()}
    as.data.frame(get_results(read_excel(path = input$file$datapath[input$file$name==input$Select], skip = 1, trim_ws = TRUE)))
    
  })
  # data =  read_excel(file=input$file$datapath[input$file$name==input$Select], skip = 1, trim_ws = TRUE)
  # get_results(data)
  # })
  
  dat <- reactive({
    if(is.null(input$file)){return()}
    data = read_excel(path = input$file$datapath[input$file$name==input$Select], skip = 1, trim_ws = TRUE)
    data = drop_na(data)
    #data[,-1] = data[,-1]-min(data[,-1])
    
    return(data)
  })
  
  output$plot <- renderPlot({ 
    data = dat()
    
    cell = unname(unlist(data[ ,input$num]))   #columns of data, user selected
    my_min = min(cell)
    cell_norm = cell - my_min
    my_max = max(cell_norm)
    cell_norm = cell_norm / my_max
    logical = as.logical(get_spike_inx(cell))
    
    
    
    ggplot(data, aes(x = as.vector(as.matrix(data[,1])), y = cell_norm))+
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
      zip(file,files)
    }
  )
  
  
})