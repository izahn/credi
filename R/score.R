#' Score CREDI response data
#'
#' Employs Maximum A Posteriori (MAP) estimation to score CREDI response data from the long form.
#' @param reverse_code (Logical) Defaults to TRUE. If TRUE, then reverse coding is automated to appropriately handle the negatively worded items LF9, LF102, LFMH1, LFMH2, LFMH3, LFMH4, LFMH5, LFMH7, LFMH8, & LFMH9. If FALSE, then no reverse coding is applied.
#' @param log_file (Character) Defaults to "logfile.txt". Name for a *.txt log file providing more details on missing data, error messages, & other information.
#' @keywords CREDI
#' @export
#' @examples
#' score()



score<-function(reverse_code = TRUE, log_file = "logfile.txt"){

    # Load required pacakges

    if(!require(stats)){install.packages("stats")}; library("stats")
    if(!require(svDialogs)){install.packages("svDialogs")}; library("svDialogs")

    # Created log file
    if(!endsWith(log_file, ".txt")){log_file = paste(log_file,".txt",sep = "")}
    sink(file = log_file)
    writeLines("")
    writeLines("------------------------------------")
    writeLines("CREDI Scoring - Log File")
    writeLines(as.character(Sys.time()))
    writeLines("------------------------------------")
    writeLines("\n")
    sink()

    # Load in the response data
    out_dlgOpen = dlgOpen(title = "Select the .csv file with the CREDI response data",
                          filters = c("csv", "*.csv"))
    csv_file = out_dlgOpen$res
    if (!endsWith(csv_file, ".csv")){stop("Selected file is not a .csv file.", call. = FALSE)}
    csv_wd = paste(strsplit(csv_file,"/")[[1]][-length(strsplit(csv_file,"/")[[1]])],collapse = "/")
    setwd(csv_wd)
    input_df = read.csv(file = csv_file)

    # Clean the input data
    list_cleaned = clean(input_df = input_df, mest_df = mest_df, csv_wd = csv_wd, reverse_code = TRUE, log_file = log_file)
    cleaned_df = list_cleaned$cleaned_df
    items_noresponse = list_cleaned$items_noresponse

    # Crate data matricies
    X = model.matrix(~1 + I( (AGE-18)/10.39 ) + I( ((AGE-18)/10.39)^2 ) + I( ((AGE-18)/10.39)^3 ), data = cleaned_df)
    Y = as.matrix(cleaned_df[,-match(c("ID","AGE",items_noresponse), names(cleaned_df))]); Y[is.na(Y)] = -9L
    MU_LF = X%*%as.matrix(B) #NxK (matrix)
    MU_SF = X%*%as.numeric(beta) #Nx1

    # Obtain necessary parameter matricies
    inds_exclude = match(items_noresponse, mest_df$CREDI_code)
    LAMBDA = as.matrix(mest_df[-inds_exclude,c("MOT","COG","LANG","SEM")])
    TAU = as.vector(mest_df$tau[-inds_exclude])
    ALPHA = as.vector(mest_df$alpha[-inds_exclude])
    DELTA = as.vector(mest_df$delta[-inds_exclude])

    # Obtain necessary constants
    J = ncol(Y);
    K = 4L
    P = 3L
    N = as.integer(nrow(Y))
    invS = as.matrix(invS)
    SIGMA_SQ= exp(X%*%as.numeric(gamma))


    # initialize the theta values
    THETA0_LF = MU_LF #NxK (matrix)
    THETA0_SF = MU_SF #Nx1 (matrix)

    # Conduct the optimization
    MAP_LF = 0.*THETA0_LF + NA
    MAP_SF = 0.*THETA0_SF + NA
    SE_LF = MAP_LF
    SE_SF = MAP_SF
    writeLines("\nScoring:")
    pb<-txtProgressBar(min = 0, max = N, initial = 0, style = 3)
    for (i in 1:N){

      # Score the long form
      out_LF = optim(par = as.vector(THETA0_LF[i,]),
                     fn = lf_posterior_density,
                     gr = lf_grad_posterior_density,
                     Yi = as.vector(Y[i,]),
                     MUi = as.vector(MU_LF[i,]),
                     invS =invS,
                     TAU = TAU,
                     LAMBDA = LAMBDA,
                     J = J,
                     K = K,
                     method = "BFGS",
                     hessian = TRUE)
      if(out_LF$convergence == 0){
        MAP_LF[i,] = out_LF$par
        fisherInfo = out_LF$hessian
        SE_LF[i,] = sqrt(diag(solve(fisherInfo,diag(K))))
      }


      # Score the short form
      out_SF = optim(par = as.vector(THETA0_SF[i,]),
                     fn = sf_posterior_density,
                     Yi = as.vector(Y[i,]),
                     MUi = as.vector(MU_SF[i,]),
                     SIGMA_SQi = as.numeric(SIGMA_SQ[i]),
                     DELTA = as.vector(DELTA),
                     ALPHA = as.vector(ALPHA),
                     J = as.integer(J),
                     method = "BFGS",
                     hessian = TRUE)
      if(out_SF$convergence == 0){
        MAP_SF[i,] = out_SF$par
        SE_SF[i,] = sqrt(1.0/out_SF$hessian)
      }

      setTxtProgressBar(pb, i)
    }

    # Clean up the MAP_LF and SE_LF
    MAP_LF = data.frame(round(MAP_LF,3))
    SE_LF = data.frame(round(SE_LF,3)); names(SE_LF) = paste(names(SE_LF),"_SE", sep = "")

    # Clean up the MAP_SF and SE_SF
    MAP_SF = data.frame(OVERALL = round(MAP_SF,3)+50)
    SE_SF = data.frame(OVERALL_SE = round(SE_SF,3))


    # Put in the input
    output_df = cbind(data.frame(ID = cleaned_df$ID), MAP_LF, MAP_SF, SE_LF, SE_SF)

    # Write out the data
    out_dlgDir = dlgSave(default = csv_wd, title = "Save as", gui = .GUI)
    out_csv = paste(strsplit(out_dlgDir$res,"/")[[1]],collapse = "/")
    if (!endsWith(out_csv,".csv")){out_csv = paste(out_csv, ".csv", sep = "")}
    write.csv(output_df, file = out_csv, row.names = FALSE)


    setwd(csv_wd)
    sink(file = log_file, append = TRUE)
    writeLines("\nScoring successful.")
    sink()
}
