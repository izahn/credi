#' Score CREDI data
#'
#' Score CREDI data and save results as a .csv file.
#' @keywords CREDI
#' @export
#' @examples
#' score()

score<-function(){

    # Load required pacakges

  if(!require(stats)){install.packages("stats")}; library("stats")
  if(!require(svDialogs)){install.packages("svDialogs")}; library("svDialogs")

  # Load in the response data
  out_dlgOpen = dlgOpen(title = "Select the .csv file with the CREDI response data",
                        filters = c("csv", "*.csv"))
  csv_file = out_dlgOpen$res
  if (!endsWith(csv_file, ".csv")){stop("Selected file is not a .csv file.")}
  csv_wd = paste(strsplit(csv_file,"/")[[1]][-length(strsplit(csv_file,"/")[[1]])],collapse = "/")
  setwd(csv_wd)
  input_df = read.csv(file = csv_file)
  N = nrow(input_df)

  # Make sure the response data is correct format
  names(input_df) = toupper(names(input_df))
  if (!unique(names(input_df)%in%input_names)){stop("Unknown variable names in csv files")}

  # Clean the ID variable as needed
  if (!"ID"%in%names(input_df)){input_df = transform(input_df, ID = 1:N)}
  if (sum(is.na(input_df$ID))>0){stop("All rows must have a unique identifier.")}
  if (length(unique(input_df$ID))<N){stop("Identifer is not unique identifier.")}

  # Clean the MONTHS as needed
  if (sum(is.na(input_df$MONTHS))>0){stop("MONTHS must be specified for each row.")}


  # Clean the missing variables
  input_df[is.na(input_df)] = -999L

  # Define needed variables
  Y = as.matrix(input_df[,-match(c("ID","MONTHS"),names(input_df))]) #NxJ (matrix)
  X = model.matrix(~-1 + I( (MONTHS-18)/10.39 ) +
                     I( ((MONTHS-18)/10.39)^2 ) + I( ((MONTHS-18)/10.39)^3 ),
                   data=input_df) #NxP (matrix)
  MU = X%*%B #NxK (matrix)


  # initialize the theta values
  THETA0 = MU #NxK (matrix)


  # Conduct the optimization for scoring
  MAP = 0.*THETA0 + NA
  MAP_SE = MAP
  pb<-txtProgressBar(min = 0, max = N, initial = 0, style = 3)
  for (i in 1:N){
    out = optim(par = as.vector(THETA0[i,]),
                fn = posterior_density,
                gr = grad_posterior_density,
                Yi = as.vector(Y[i,]),
                MUi = as.vector(MU[i,]),
                invS =invS,
                TAU = TAU,
                LAMBDA = LAMBDA,
                J = J,
                K = K,
                method = "BFGS",
                hessian = TRUE)
    if(out$convergence == 0){
      MAP[i,] = out$par

      fisherInfo = out$hessian
      MAP_SE[i,] = sqrt(diag(solve(fisherInfo,diag(K))))
    }
    setTxtProgressBar(pb, i)
  }

  # Clean up the MAP and MAP_SE
  MAP = data.frame(round(MAP,3))
  MAP_SE = data.frame(round(MAP_SE,3)); names(MAP_SE) = paste(names(MAP_SE),"_SE", sep = "")

  # Put in the input
  output_df = cbind(data.frame(ID = input_df$ID), MAP, MAP_SE)

  # Write out the data
  out_dlgDir = dlgSave(default = csv_wd, title = "Save as", gui = .GUI)
  out_csv = paste(strsplit(out_dlgDir$res,"/")[[1]],collapse = "/")
  if (!endsWith(out_csv,".csv")){out_csv = paste(out_csv, ".csv", sep = "")}
  write.csv(output_df, file = out_csv, row.names = FALSE)

  return(output_df)

}
