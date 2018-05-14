#' Clean input CREDI response data
#'
#' This calculates the posterior density function.
#' @param input_df Raw input data from the user-specified csv file.
#' @param mest_df Measurement parameter estimates.
#' @param csv_wd Working directory where input *.csv file is located.
#' @param reverse_code Logical. If TRUE, then reverse codes LF9, LF102, LFMH1, LFMH2, LFMH3, LFMH4, LFMH5, LFMH7, LFMH8, & LFMH9. If FALSE, then no reverse coding is applied.
#' @param log_file Name of the *.txt log file for the CREDI scoring.
#' @keywords CREDI
#' @export
#' @examples
#' clean(input_df, mest_df, csv_wd, reverse_code,log_file)      


 clean<-function(input_df, mest_df, csv_wd, reverse_code,log_file){    
      # Input: 
      #  input_df - User defined input file, with: 
      #                 a) LF item naming convention needed.
      #                 b) Unique identifier needed.
      #                 c) Age in month needed
      #  mest_df - Item parameters and other information like naming conversion and reverse coding.
      #  csv_wd - File location with the data
      #  reverse_code - (logical) If TRUE then reverse codes long form items that are negatively worded. 
      #                           Otherwise does not implement reverse coding.
      # Output: List with the following:
      #  cleaned_df - Cleaned data (i.e., missing data codes, reverse coding as necessary)
      #  items_noresponse - character vector with items missing all responses
      
      setwd(csv_wd)
      

      # Ensure the naming of the response data is in the correct format
      input_names = c("ID","AGE",mest_df$CREDI_code)
      names(input_df) = toupper(names(input_df))
      cols_vnames = which(!names(input_df)%in%input_names)
      if (length(cols_vnames)>0){
        stop_message = paste("Unknown variable names in csv file: ", paste(names(input_df)[cols_vnames], collapse = ", "), sep = "")
        sink(file = log_file, append = TRUE)
        writeLines(stop_message)
        sink()
        stop(stop_message)
      }
      
      # Ensure that there is a unique ID variable for each observations
      if (!"ID"%in%names(input_df)){
        stop_message = "*Error: An identifier variable named ID must be included in the .csv file."
        sink(file = log_file, append = TRUE)
        writeLines(stop_message)
        sink()
        stop(stop_message, call. = FALSE)
      }
      if (sum(is.na(input_df$ID))>0){
        stop_message = "*Error: ID variable missing for some observations."
        sink(file = log_file, append = TRUE)
        writeLines(stop_message)
        sink()
        stop(stop_message, call. = FALSE)
      }
      if (length(unique(input_df$ID))<nrow(input_df)){
        stop_message = "*Error: Identifer is not unique identifier. See log file more more details."
        sink(file = log_file, append = TRUE)
        writeLines(stop_message)
        writeLines("Duplicate identifiers include:")
        writeLines( paste("ID = ", paste(input_df$ID[which(duplicated(input_df$ID))], collapse = ", "), sep = "" ) )
        sink()
        stop(stop_message, call. = FALSE)
      }
      
      # Check that all variables outside of ID are numeric
      not_numeric = NULL
      for (j in seq(1,ncol(input_df))){
        if( !is.numeric(input_df[,j]) & names(input_df)[j]!="ID" ){
          not_numeric = c(not_numeric,names(input_df)[j])
        }
      }
      if (!is.null(not_numeric)){
        stop_message = "*Error: AGE and all item response variables must be in numeric format. At least one of these variables may contain non-numeric values. See log file more more details."
        sink(file = log_file, append = TRUE)
        writeLines(stop_message)
        writeLines("Check that the following variable(s) contain only numeric values in the .csv file:")
        writeLines( paste(not_numeric, collapse = ", ") )
        sink()
        stop(stop_message, call. = FALSE)
      }
      
      
      # Create a log of number of observations that must be discarded
      N_input = nrow(input_df)
      discard_df = data.frame(Reason = rep(NA,3), Number = rep(NA,3))
      dr = 0
      
      # Clean the AGE as needed
      #input_df$AGE[c(2, 16, 28, 45:50)] = NA # works
        # Missing age
        rows_mi_age = which(is.na(input_df$AGE))
        if (length(rows_mi_age)>0){
          sink(file = log_file, append = TRUE)
          writeLines(paste("\n* Warning:\n  The following ", length(rows_mi_age)  , " observation(s) are missing AGE values and cannot be scored:\n  ID = ", paste(input_df$ID[rows_mi_age], collapse = ", "), sep = ""))
          sink()
          dr = dr+1; discard_df$Reason[dr] = "Missing age values"; discard_df$Number[dr] = length(rows_mi_age)
          input_df = input_df[-rows_mi_age, ]
        }
        
        # Age outside of range.
        if (nrow(input_df)>0){
          rows_out_age = which(input_df$AGE<0 | input_df$AGE>36)
          if (length(rows_out_age)>0){
            sink(file = log_file, append = TRUE)
            writeLines(paste("\n* Warning:\n  The following ", length(rows_out_age) ," observation(s) cannot be scored because the AGE values are outside of 0-36 months:\nID = ", paste(input_df$ID[rows_out_age], collapse = ", "), sep = ""))
            sink()
            dr = dr+1; discard_df$Reason[dr] = "Age values outside of 0-36 months"; discard_df$Number[dr] = length(rows_out_age)
            input_df = input_df[-rows_out_age, ]
          }
        }
      
        if (nrow(input_df)>0){
          # Clean the missing responses data
          cols_Q = which(startsWith(names(input_df), "LF"))
          temp_df = input_df[,cols_Q]
          temp_df[temp_df!=0L & temp_df!=1L] = NA
        
          # Discard observations with fewer than 5 item responses
          num_nonmi_y = apply(temp_df, MARGIN = 1L, function(X){sum(!is.na(X))})
          rows_toofew_y = which(num_nonmi_y<5L)
          if (length(rows_toofew_y)>0){
            sink(file = log_file, append = TRUE)
            writeLines(paste("\n* Warning:\n  The following ", length(rows_toofew_y) ," observation(s) contain less than 5 non-missing item responses and will not be scored:\n  ID = ", paste(input_df$ID[rows_toofew_y], collapse = ", "), sep = ""))
            sink()
            dr = dr+1; discard_df$Reason[dr] = "Less than 5 item responses"; discard_df$Number[dr] = length(rows_toofew_y)
            input_df = input_df[-rows_toofew_y, ]
          }
        }
        
        # Check if there are any remaining observations to score
        if (nrow(input_df)==0){
          stop_message = paste("\n* Error:\n  All ", N_input," observations have been discarded for the following reason(s): \n", sep = "")
          sink(file = log_file, append = TRUE)
          writeLines(stop_message)
          print(discard_df[complete.cases(discard_df), ])
          writeLines("* Error: No observations to score.")
          sink()
          writeLines(stop_message)
          print(discard_df[complete.cases(discard_df), ])
          stop("No observations to score.", call. = FALSE)
        }
        
        # Check if they want to continue given the discarding
        if(nrow(input_df)<N_input){
          discard_df = discard_df[complete.cases(discard_df), ]
          discard_df = transform(discard_df, Percent = paste(round(100*Number/N_input,1),"%",sep = ""))
          N_discarded = sum(discard_df$Number); Pct_discarded = round(100*N_discarded/N_input,1)
          
          sink(file = log_file, append = TRUE)
          writeLines(paste("\n* Warning:\n  A total of ", N_discarded, " (", Pct_discarded,"%) observation(s) cannot be scored for the following reason(s):",sep =""))
          print(discard_df[complete.cases(discard_df), ])
          sink()
          
          writeLines(paste("\n* Warning:\n  A total of ", N_discarded, " (", Pct_discarded,"%) observation(s) cannot be scored (see logfile for more details):",sep =""))
          print(discard_df[complete.cases(discard_df), ])
          

          x<-as.character(readline(prompt = "Would you like to continue? [Y/N]: "))
          x <- toupper(x)
          cut = 0
          while(cut == 0){
            
            if (x == "Y"){
              cut = 1;
            } else if (x == "N"){ 
              cut = 1
              stop("Scoring canceled.", call. = FALSE) 
            } else {
              x<-as.character(readline(prompt = "Would you like to continue? [Y/N]:"))
              x <- toupper(x)
              cut = 0
            }
            
          } #end while
        
        } # end if
        
        
        # Create the cleaned_df version 
        N = nrow(input_df)
        cleaned_df = data.frame(mat.or.vec(nr = N, nc = nrow(mest_df)+2)+NA)
        names(cleaned_df) = c("ID","AGE", mest_df$CREDI_code)
        cleaned_df$ID = input_df$ID; cleaned_df$AGE = input_df$AGE
        for (j in cols_Q){
          
          # Only bring in answers that are 1 or 0. Otherwise keep missing
          col_j = which( names(cleaned_df) %in% names(input_df)[j] )
          rows_j = which(input_df[,j] == 1L | input_df[,j]==0)
          cleaned_df[rows_j,col_j] = input_df[rows_j,j]
          
        }
      
        
        # Print out missing data descriptive statistics
        miss_df = data.frame(round(100*apply(cleaned_df[,-c(1,2)], 2, function(X){sum(is.na(X))})/N,2))
        names(miss_df) = c("Pct_Missing")
        miss_df = subset(miss_df, Pct_Missing>0)
        if (nrow(miss_df)>0){
          
          items_noresponse = row.names(miss_df)[miss_df$Pct_Missing==100]
          
          inds_order = sort(miss_df$Pct_Missing, decreasing = TRUE, index.return = TRUE)
          miss_df2 = data.frame(Item = row.names(miss_df)[inds_order$ix], Pct_Missing = miss_df$Pct_Missing[inds_order$ix])
          
          writeLines(paste("\nThe following items contained missing responses from all individuals: ", 
                           paste(subset(miss_df2, Pct_Missing==100)$Item, collapse = ", "), sep = ""))
          
          miss_df3 = subset(miss_df2, Pct_Missing<100)
          miss_df3 = transform(miss_df3, Pct_Missing = round(Pct_Missing, 1))
          miss_df3 = transform(miss_df3, Pct_Missing = paste(Pct_Missing,"%", sep = ""));
          
          #miss_df4 = data.frame(mat.or.vec(nr = 1, nc = nrow(miss_df3))+NA)
          #names(miss_df4) = miss_df3$Item
          #miss_df4[1,] = miss_df3$Pct_Missing
          sink(file = log_file, append = TRUE)
          writeLines("\nMissingness rates of items responses:")
          print(miss_df2)
          sink()
          
          
        }
        
        if (reverse_code == TRUE){
            # Reverse code as needed
            items_reverse = mest_df$CREDI_code[mest_df$RevCoded]
            cols_reverse = which(!is.na(match(names(cleaned_df),items_reverse)))
            for (j in cols_reverse){
              rows_j = which(cleaned_df[,j]==1L | cleaned_df[,j]==0L)
              cleaned_df[rows_j,j] = as.integer(1L-cleaned_df[rows_j,j])
            }
        }
        
        out_list = list(cleaned_df = cleaned_df, 
                        items_noresponse = items_noresponse)
        
        return(out_list)
   
 }     
        