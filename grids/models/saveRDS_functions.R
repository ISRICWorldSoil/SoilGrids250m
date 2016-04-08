## based on: http://stackoverflow.com/questions/28927750/what-is-the-method-to-save-objects-in-a-compressed-form-using-multiple-threads-c

## Save RDS in parallel
saveRDS.xz <- function(object,file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pxz -T",threads," > ",file),"wb")
  saveRDS(object, file = con, compress=FALSE)
  close(con)
}

readRDS.xz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pxz -d -k -c -T",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

saveRDS.gz <- function(object,file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

readRDS.gz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

readRDS.p <- function(file,threads=parallel::detectCores()) {
  #Hypothetically we could use initial bytes to determine file format, but here we use the Linux command file because the readBin implementation was not immediately obvious
  fileDetails <- system2("file",args=file,stdout=TRUE)
  selector <- sapply(c("gzip","XZ"),function (x) {grepl(x,fileDetails)})
  format <- names(selector)[selector]
  if (format == "gz") {
    object <- readRDS.gz(file, threads=threads)
  } else if (format == "XZ") {
    object <- readRDS.xz(file, threads=threads)
  } else {
    object <- readRDS(file)
  }
  return(object)
}