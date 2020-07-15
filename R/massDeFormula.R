#'  molecular mass from chemical formula 
#'
#'  Calculate molecular mass based on atomic composition    
#' 
#' @param comp (character) atomic compostion 
#' @param massTy (character) 'mono' or 'average'
#' @param rmEmpty (logical) suppress empty entries
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of messages produced
#' @return numeric vector with mass
#' @seealso \code{\link[wrMisc]{convToNum}}
#' @examples
#' massDeFormula(c("12H12O","HO"," 2H 1 Se, 6C 2N","HSeCN"," ","e"))
#' @export
massDeFormula <- function(comp,massTy="mono",rmEmpty=FALSE,silent=FALSE,callFrom=NULL){
  ## calculate molecular mass based on composition formula (sum formula: number & element)
  ## 'comp' .. character vector with molecular composition(s)
  ##
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="massDeFormula")
  if(length(wrMisc::naOmit(comp)) < length(comp)){
    if(!silent) message(fxNa,sum(is.na(comp))," entries of 'comp' are NA  (remove) !")
    comp <- comp[wrMisc::naOmit(match(wrMisc::naOmit(comp), comp))] }
  ## rm NAsclean heading space, 
  comp <- gsub("^ +","",gsub(" +$","",comp))
  if(rmEmpty) {
    if(any(comp=="")) {
      if(!silent) message(fxNa,sum(comp=="")," some entries of 'comp' are empty  (remove) !")
      comp <- comp[which(comp !="")] }    
  }
  msg <- "Can't find any element names (must start with caps,'e' or 'z')'"
  chEm <- nchar(comp) <1
  if(any(chEm)) comp[which(chEm)] <- "z" 
  El <- up <- gregexpr("[[:upper:]]|e|z",comp)  
  chMaj <- sapply(up, function(x) any(x <1))
  if(any(chMaj)) stop(msg," in ",comp[which(chMaj)])
  ay <- gregexpr("[[:lower:]]", comp)
  for(i in which(sapply(ay, function(x) any(x >0)))) {            # correct if upper followed by lower caps
    tmp <- which(El[[i]] %in% (ay[[i]]-1))
    El[[i]][tmp] <- El[[i]][tmp] +1 }
  form <- list()
  for(i in 1:length(El)) {
    begStr <- c(1,El[[i]][-length(El[[i]])] +1)           # beginning of string with number & element
    y <- substring(comp[[i]], begStr,El[[i]])              # isolated series until capital/lower letter
    sig <- grep("-",y)
    num <- gsub("[[:alpha:]]|[[:blank:]]|[[:punct:]]","",y)
    num[which(num=="")] <- "1"
    num <- as.numeric(num)
    if(length(sig) >0) num[sig] <- -1*num[sig]
    if(length(sig) <length(y)) num[-1*sig] <- paste("+",num[-1*sig],sep="")  # add '+'
    form[[i]] <- matrix(c(num,substring(comp[[i]],up[[i]],El[[i]])), ncol=2, dimnames=list(NULL,c("n","elem"))) }
  names(form) <- comp
  ## convert extracted/cleaned sum formula in mass :
  atMa <- .atomicMasses()[,massTy=massTy]
  tmp <- lapply(form,function(x) x[,2] %in% names(atMa))
  chEl <- sapply(tmp,sum)
  usePep <- 1:length(comp)
  if(any(chEl < sapply(form,nrow))) {
    if(all(chEl <1)) stop(" Can't identify any of the names given via .atomicMasses()")
    nonIdEl <- unlist(sapply(form,function(x) x[which(is.na(match(x[,2],names(atMa)))),2]))
    if(any(chEl <1) & !silent) message(fxNa, " can't find ",wrMisc::pasteC(sub("z"," ",nonIdEl),quote="'")," .. setting to 0 mass")
    corEl <- which(chEl < sapply(form,nrow)) 
    form[corEl] <- lapply(form[corEl], function(x) {x[which(x[,2] %in% nonIdEl),1] <- "0"; x})
  }
  mass <- sapply(form,function(x) sum(as.numeric(x[,1])*atMa[match(x[,2],names(atMa))]))
  names(mass) <- sapply(form,function(x) paste(paste(x[,1],x[,2],sep=""),collapse=""))
  chNa <- names(mass) =="0z"
  if(any(chNa)) {y <- which(chNa); mass[y] <- 0; names(mass)[y] <- ""}
  mass }
  
 #' @export
.atomicMasses <- function() {
  ## return matrix of atomic masses : 1st col for average mass and 2nd col for mono-isotopic
  ## based on http://www.ionsource.com/Card/Mass/mass.htm in ref to http://physics.nist.gov/Comp   (~agree in http://www.weddslist.com/ms/tables.html)
  mass <- cbind(aver=c(1.007940, 12.010700, 14.006700, 15.999400, 30.973761, 32.065000, 78.960000, 5.48579909e-4,0),
    mono=c(1.0078250321, 12, 14.0030740052, 15.9949146221, 30.97376151, 31.97207069, 79.9165196, 5.48579909e-4,0))
  rownames(mass) <- c("H","C","N","O","P","S","Se","e","ze")
  mass }
   
