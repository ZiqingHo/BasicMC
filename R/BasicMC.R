#'@title (General)Convert raw data into a transition matrix
#'@name createtrMat
#'
#'@description General function to create a transition matrix from raw data
#'
#'@param y Data frame that contains the frequency of the transitions of each state/Raw data in the form of numeric vector, character vector or matrix.
#'@param states The character or numeric vector of states.
#'@param transition Choose whether the data uploaded is in the form of transition matrix or not. The default value is FALSE.
#'
#'@return The transition matrix of the Markov chain
#'
#'@examples
#'#raw data in "data.frame" class is used
#'data(bigram)
#'statenames<-c("a", "b", "c", "d", "e", "f", "g", "h", "i",
#'              "j", "k", "l", "m", "n", "o", "p", "q", "r", "s",
#'              "t", "u", "v", "w", "x", "y", "z")
#'createtrMat(bigram,statenames)
#'
#'#vector of raw data in "numeric" class is used
#'x <- c(1,4,1,4,3,4,4,4,3,1,3,2,3,3,3,4,2,2,3)
#'createtrMat(x,c(1:4))
#'
#'#vector of raw data in "character" class is used
#'y<-c("a","d","a","d","c","d","d","d","c",
#'     "a","c","b","c","c","c","d","b","b","c")
#'createtrMat(y,c("a","b","c","d"))
#'
#'#data in "matrix" class is used
#'y<-matrix(c(5,6,7,8),byrow=TRUE,nrow=2)
#'createtrMat(y,c(1,2))
#'
#'@export
createtrMat<-function(y,states,transition=FALSE){
  if(transition==FALSE){
  if(is.vector(y) == TRUE){
    len<-length(y)
    trMat<-matrix(nrow = length(states), ncol = length(states), 0)

    if(is.numeric(y) == TRUE){
      for(i in 1:(len - 1)){
        trMat[y[i],y[i+1]]<-trMat[y[i],y[i+1]] + 1
      }
      for(j in 1:length(states)){
        total<-sum(trMat[j,])
        trMat[j,]<-trMat[j,]/total
      }
    }

    if(is.character(y) == TRUE){
      l1<-c()
      for(i in 1:length(y)){
        for(j in 1:length(states)){
          if(y[i] == states[j]){
            state<-j
            l1<-append(l1,state)
          }
        }
      }
      for(i in 1:(len - 1)){
        trMat[l1[i],l1[i+1]]<-trMat[l1[i],l1[i+1]] + 1
      }
      for(j in 1:length(states)){
        total<-sum(trMat[j,])
        trMat[j,]<-trMat[j,]/total
      }
    }
    nr1<-nrow(trMat)
  }


  if(is.data.frame(y) == TRUE){
    nr<-nrow(y)
    nc<-ncol(y)
    trMat <- matrix(c(rep(c(rep(0,nr)),nr)),byrow=TRUE,nrow=nr)

    nr1<-nrow(trMat)

    for(i in 1:nr1){
      total<-sum(y[i,-1])
      for(j in 2:nc){
        trMat[i,j-1]<-y[i,j]/total
      }
    }
  }

  if(is.matrix(y) == TRUE){
    nr<-nrow(y)
    nc<-ncol(y)
    trMat <- matrix(c(rep(c(rep(0,nr)),nr)),byrow=TRUE,nrow=nr)

    nr1<-nrow(trMat)

    for(i in 1:nr1){
      total<-sum(y[i,])
      for(j in 1:nc){
        trMat[i,j]<-y[i,j]/total
      }
    }
  }

  if(length(states) != nr1){
    stop("The length of argument states must be equal to the number of states.")
  }

  colnames(trMat)<-c(states)
  rownames(trMat)<-c(states)

  class(trMat)<-"MarkovChain"
  return(trMat)
  }else{

    if(is.matrix(y) == TRUE){
      nr<-nrow(y)
      nc<-ncol(y)
      trMat <- matrix(c(rep(c(rep(0,nr)),nr)),byrow=TRUE,nrow=nr)

      for(i in 1:nr){
        for(j in 1:nc){
          trMat[i,j]<-y[i,j]
        }
      }
    }

    if(is.data.frame(y) == TRUE){
      nr<-nrow(y)
      nc<-ncol(y)
      trMat <- matrix(c(rep(c(rep(0,nr)),nr)),byrow=TRUE,nrow=nr)


      for(i in 1:nr){
        for(j in 2:nc){
          trMat[i,j-1]<-y[i,j]
        }
      }
    }
    if(length(states) != nrow(trMat)){
      stop("The length of argument states must be equal to the number of states.")
    }
    colnames(trMat)<-c(states)
    rownames(trMat)<-c(states)
    class(trMat)<-"MarkovChain"
    return(trMat)
  }

}
#'@title Convert the data into a transition matrix for stock market
#'
#'@name stock
#'
#'@description This method is used to create a transition matrix from raw data of stock prices.
#'
#'@param y Data of the stock prices
#'@param pricecolumn The index of the column containg the prices
#'@param interval The range for the state "non-significant". The default value is 10. If interval = 0, "non-significant" state will exist when no price change.
#'
#'@return The transition matrix of the Markov chain
#'@examples
#'data(ftse)
#'stock(ftse,3,5)

#'@export

stock<-function(y,pricecolumn,interval=10){

  if(is.numeric(pricecolumn)==FALSE){
    stop("The argument pricecolumn must be the index of price column.")
  }

  nr<-nrow(y)
  statelist<-c()

  for(i in 2:nr){
    pricechange<-y[i,pricecolumn]-y[i-1,pricecolumn]

    if(pricechange>interval){
      state<-1
    }else if(pricechange<(-interval)){
      state<-3
    }else if(pricechange<interval|pricechange>(-interval)){
      state<-2
    }

    statelist<-append(statelist,state)
  }
  n<-length(statelist)
  trMat <- matrix(c(rep(c(0,0,0),3)),byrow=TRUE,nrow=3)

  nr1<-nrow(trMat)
  nc<-ncol(trMat)

  for(j in 2:n){
    for(k in 1:3){
      if(statelist[j-1]==k){
        trMat[k,statelist[j]]<-trMat[k,statelist[j]]+1
      }
    }
  }

  for(l in 1:3){
    total<-sum(trMat[l,])
    trMat[l,1]<- trMat[l,1]/total
    trMat[l,2]<- trMat[l,2]/total
    trMat[l,3]<- trMat[l,3]/total
  }

  colnames(trMat)<-c("Rise","Non-significant","Fall")
  rownames(trMat)<-c("Rise","Non-significant","Fall")


  if(is.nan(trMat[2,1]) == TRUE){
    trMat<-trMat[-2,-2]
  }

  class(trMat)<-"MarkovChain"
  return(trMat)

}

#'@title Convert the data into a transition matrix for weather forecasting
#'
#'@name weather
#'
#'@description This method is used to create a transition matrix from the data of weather
#'
#'@param y Data of the weather
#'@param weathercolumn The index of the column containing weather types
#'@param sunny The character vector for 'sunny' state
#'@param cloudy The character vector for 'cloudy' state
#'@param rainy The character vector for 'rainy' state
#'@param other The character vector for other weather types. The default value = NULL.
#'
#'@return The transition matrix of the Markov chain
#'
#'@examples
#'data(wsoton)
#'rainnystate<-c("Light rain", "Moderate rain", "Patchy rain possible",
#'              "Light drizzle", "Light rain shower")
#'weather(wsoton,2,c("Sunny"),c("Cloudy","Partly cloudy"),rainnystate,"Overcast")

#'@export
weather<-function(y,weathercolumn,sunny,cloudy,rainy,other=NULL){

  if(is.numeric(weathercolumn)==FALSE){
    stop("The argument weathercolumn must be the index of weather column.")
  }

  nr<-nrow(y)
  statelist<-c()
  state<-NULL

  for(i in 1:nr){

    for(j in 1:length(sunny)){
      if(y[i,weathercolumn]==sunny[j]){
        state<-1
      }
    }
    for(k in 1:length(cloudy)){
      if(y[i,weathercolumn]==cloudy[k]){
        state<-2
      }
    }
    for(l in 1:length(rainy)){
      if(y[i,weathercolumn]==rainy[l]){
        state<-3
      }
    }

    if(!is.null(other)){
      for (m in 1:length(other)){
        if(y[i,weathercolumn]==other[m]){
          state<-m+3
        }
      }
    }
    statelist<-append(statelist,state)
  }

  n<-length(statelist)
  n1<-3+length(other)
  trMat <- matrix(c(rep(c(rep(0,n1)),n1)),byrow=TRUE,nrow=n1)

  nc<-ncol(trMat)

  for(j in 2:n){
    for(k in 1:nc){
      if(statelist[j-1]==k){
        trMat[k,statelist[j]]<-trMat[k,statelist[j]]+1
      }
    }
  }

  for(l in 1:nc){
    total<-sum(trMat[l,])
    for(m in 1:nc){
      trMat[l,m]<- trMat[l,m]/total}

  }
  colnames(trMat)<-c("Sunny","Cloudy","Rainy",c(other))
  rownames(trMat)<-c("Sunny","Cloudy","Rainy",c(other))

  list1<-c()
  for(i in 1:nc){
    if(is.nan(trMat[i,i]) == TRUE){
      list1<-append(list1,i)
    }else{
      list1
    }
  }

  trMat<-trMat[-list1,-list1]
  class(trMat)<-"MarkovChain"
  return(trMat)

}

#`Set validity and initialize transition matirx`
#'@import methods
#'@title Check and initialise the transition matrix used
#'
#'@name propertrMat
#'
#'@description This method is used to check and initialise the transition matrix used
#'
#'@param trMatrix Transition matrix of the Markov chain
#'
#'@return Transition matrix of the Markov chain
#'@references Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
#'@examples
#'egstate <- c("a", "b", "c")
#'MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#'propertrMat(MC)

#'@export
propertrMat <- function(trMatrix){

  for(i in 1:ncol(trMatrix)){
  if(sum(trMatrix[i,])<0.9999999999){
    stop("The sum of transition probabilities from a state to all other states must be 1.")}
}
  if(!all(trMatrix>=0)|!all(trMatrix<=1)){
    stop("The transition probability must be in the interval [0,1].")}

  if(ncol(trMatrix)!=nrow(trMatrix)){
    stop("The argument trMatrix must be a square matrix.")}

  if(all(is.null(colnames(trMatrix)),is.null(rownames(trMatrix)))){
    nc <- ncol(trMatrix)
    colnames(trMatrix) <- as.character(seq_len(nc))
    rownames(trMatrix) <- as.character(seq_len(nc))
  }

  if (is.null(colnames(trMatrix))){
    colnames(trMatrix)<-rownames(trMatrix)
  }
  if(is.null(rownames(trMatrix))){
    rownames(trMatrix)<-colnames(trMatrix)
  }

  if(any(colnames(trMatrix) != rownames(trMatrix))){
    rownames(trMatrix)<-colnames(trMatrix)
  }

  if(is.data.frame(trMatrix) == TRUE){
    trMatrix<-data.matrix(trMatrix)
  }
  return(trMatrix)

}

#`Display transition matirx`
#'@title Display the transition matrix
#'
#'@name displaytrM
#'
#'@description This method is used to display the transition matrix of a finite-state discrete time-homogeneous Markov Chain.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'@rdname displaytrM
#'@return The transition matrix of the corresponding Markov chain
#'@references Chambers, J. (2008). Software for Data Analysis: Programming with R. New York: Springer Science+Business Media, LLC.
#'@examples
#'egstate <- c("a", "b", "c")
#' MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#' displaytrM(MC)
#'
#'
#' @export
displaytrM<-function(trMatrix){
  trMatrix<-propertrMat(trMatrix)
  return(trMatrix)
}

#`Display the states
#'@name displayState
#'
#'@title Display the states of a Markov chain
#'
#'@description This method is used to display the states of a transition matrix of a finite-state discrete time-homogeneous Markov Chain.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'
#'@return The character vector of the states of the corresponding Markov chain.
#'
#'@references Chambers, J. (2008). Software for Data Analysis: Programming with R. New York: Springer Science+Business Media, LLC.
#'
#'@examples
#' egstate <- c("a", "b", "c")
#' MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#' displayState(MC)
#'
#'@rdname displayState
#'
#'
#'@export
displayState<-function(trMatrix){
  trMatrix<-propertrMat(trMatrix)
  return(colnames(trMatrix))
}

#`Obtain the stationary distribution
#'@name stationary
#'
#'@title Function to obtain the stationary distribution of a Markov Chain
#'
#'@description This method is used to get the stationary distribution of a finite discrete time-homogeneous Markov chain.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'
#'@return A vector of stationary distribution corresponding to the Markov chain.
#'
#'@references Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
#'
#'@examples
#' egstate <- c("a", "b", "c")
#' MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#' stationary(MC)
#'
#'@rdname stationary
#'
#'@export
stationary<-function(trMatrix) {
  trMatrix<-propertrMat(trMatrix)

  eigenv<-eigen(t(trMatrix))
  real <- Re(eigenv$vectors[,1])
  prob<- real/sum(real)

  if(is.null(prob)){
    warning("No steady distribution exists.")
    print(NULL)
  }else{
    result<-rbind(c( ),prob)
    colnames(result)<-c(colnames(trMatrix))

    return(result)
  }
}

#`Calculate the n-step transition matrix
#'@name trProb
#'
#'@title Calculate the n-step transition matrix
#'
#'@description Function to calculate the n-step transition matrix of a finite-state discrete time-homogeneous Markov Chain.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'@param steps The number of steps
#'
#'@return A matrix of n-step transition probabilities of the states
#'@references Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
#'
#'@examples
#'egstate <- c("a", "b", "c")
#'MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#'n <- 5
#'trProb(MC,n)
#'
#'@rdname trProb
#'
#'@export
trProb<-function(trMatrix,steps){
  trMatrix<-propertrMat(trMatrix)

  if(steps%%1 !=0 |steps <=0){
    stop("The argument steps must be a positive integer.")}

  result<-trMatrix

  if(steps==1){
    class(trMatrix)<-"MarkovChain"
    return(trMatrix)
  }else{
    for(i in seq(steps-1)){
      result <- result%*%trMatrix
    }

    output <- result
    class(output)<-"MarkovChain"
    return(output)}
}

#`Calculate the conditional distribution of states given a distribution of initial states
#'@name condDist
#'
#'@title Obtain the conditional distribution of the states
#'
#'@description Function to calculate the conditional distribution of states given a distribution of initial states.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'@param initialstate A vector of distribution of the initial states.
#'@param steps The number of steps
#'
#'@return A vector of n-step transition probabilities of the states given a distribution of initial states
#'
#'@references Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
#'
#'@examples
#'egstate <- c("a", "b", "c")
#'MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#'initial <- c(0.5,0.5,0)
#'n <- 5
#'condDist(MC,initial,n)
#'
#'@rdname condDist
#'
#'@export
condDist<-function(trMatrix,initialstate=NULL,steps){
  trMatrix<-propertrMat(trMatrix)

  if(steps%%1 !=0 |steps <=0){
    stop("The argument steps must be a positive integer.")}

  if(!all(initialstate<=1)|!all(initialstate>=0)){
    stop("The probability of each state must in the interval[0,1]")
  }

  if(is.null(initialstate)){
    nc<-ncol(trMatrix)
    initialstate<-c(1,rep(0,nc-1))
  }

  if(sum(initialstate)!= 1){
    stop("The sum of probability of the initial states must equal to 1.")
  }

  result<-trMatrix

  if(steps==1){
    output <- initialstate %*% result
    return(output)
  }else{
    for(i in seq(steps-1)){
      result <- result%*%trMatrix
    }

    output <- initialstate %*% result
    return(output)}
}

#`Simulate a discrete time-homogeneous Markov chain
#'@name simMC
#'
#'@title Simulation of a finite-state discrete time-homogeneous Markov Chain
#'
#'@description This method is used to simulate a finite-state discrete time-homogeneous Markov Chain.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'@param InitStateIndex The index of the initial state in the vector of states
#'@param steps The number of iterations
#'@param presentation Decide whether the states simulated is in number or state name.The default value is TRUE.
#'
#'@return A sequence of states
#'
#'@references Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
#'
#'@examples
#'egstate <- c("1", "2", "3")
#'MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#'index <- 1
#'n <- 50
#'simMC(MC,index,n)
#'
#'@rdname simMC
#'@importFrom stats rmultinom
#'@export
simMC<-function(trMatrix,InitStateIndex=NULL,steps,presentation=TRUE){
  trMatrix<-propertrMat(trMatrix)

  if(steps%%1 !=0 |steps <=0){
    stop("The argument steps must be a positive integer.")}

  state <- numeric(steps)

  if(is.null(InitStateIndex)){
    state[1] <- 1

  }else{
    if(InitStateIndex > length(colnames(trMatrix))|InitStateIndex <= 0|InitStateIndex%%1 !=0){
      stop("The argument InitStateIndex must be a valid index in the vector of state")
    }
    state[1] <- InitStateIndex

  }

  for (i in 2:steps){
    trprob <- trMatrix[state[i-1], ]
    state[i] <- which(rmultinom(1,1,trprob) == 1)}

  if(presentation == TRUE){
    state <- state
    statelist<-colnames(trMatrix)
    return(statelist[state])
  }else{
    state <- state
    return(state)}
}

#`Obtain the transient states
#'@name transient
#'
#'@title Function to identify transient states
#'
#'@description This method is used to identify the transient states of a finite-state discrete time-homogeneous Markov Chain.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'
#'@return The transient states of a finite-state discrete time-homogeneous Markov Chain.
#'
#'@references Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
#'
#'@examples
#'egstate <- c("a", "b", "c")
#'MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#'transient(MC)
#'
#'@rdname transient
#'
#'@export
transient<-function(trMatrix) {
  trMatrix<-propertrMat(trMatrix)

  nc<-ncol(trMatrix)
  index<-c()

  result<-trMatrix
  output<-trMatrix


  for(j in seq(nc)){

    for(i in 2:50){
      result<-result%*%trMatrix
      output<- output + result}

    if(output[j,j] < 1) {
      index<-union(index,c(j))
    }else{
      index<-index
    }
  }


  if(is.null(index)){
    cat("There is no transient state exists in this Markov chain.","\n")

  }else{

    reState<-(colnames(trMatrix)[index])

    cat("Transient state(s) : ")
    cat(reState,sep=", ")}
}

#`Obtain the recurrent states
#'@name recurrent
#'
#'@title Function to identify the recurrent states
#'
#'@description This method is used to identify the recurrent states of a finite-state discrete time-homogeneous Markov Chain.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'
#'@return The recurrent states of the finite-state discrete time-homogeneous Markov Chain
#'
#'@references Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
#'
#'@examples
#'egstate <- c("a", "b", "c")
#'MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#'recurrent(MC)
#'
#'@rdname recurrent
#'
#'@export
recurrent<-function(trMatrix) {
  trMatrix<-propertrMat(trMatrix)

  nc<-ncol(trMatrix)
  index<-c()

  result<-trMatrix
  output<-trMatrix

  for(j in seq(nc)){

    for(i in 2:50){
      result<-result%*%trMatrix
      output<- output + result}

    if(output[j,j] < 1) {
      index<-index
    }else{
      index<-union(index,c(j))
    }
  }

  if(is.null(index)){
    print("There is no recurrent state exists in this Markov chain.")
  }else{

    reState<-(colnames(trMatrix)[index])

    cat("Recurrent state(s) : ")
    cat(reState,sep=", ")
    cat("\n")}
}

#`Obtain the absorbing states
#'@name absorbing
#'
#'@title Function to identify the absorbing states
#'
#'@description This method is used to identify the absorbing states of a discrete time-homogeneous Markov chain.
#'
#'@param trMatrix Transition matrix of the Markov chain
#'
#'@return The absorbing states of the finite-state discrete time-homogeneous Markov Chain
#'
#'@references Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
#'
#'@examples
#'egstate <- c("a", "b", "c")
#'MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'          nrow = 3, byrow = TRUE,
#'          dimnames = list(egstate, egstate))
#'
#'absorbing(MC)
#'
#'@rdname absorbing
#'
#'@export
absorbing<-function(trMatrix) {
  trMatrix<-propertrMat(trMatrix)

  nc<-ncol(trMatrix)
  index<-c()

  for(i in 1:nc){
    if(trMatrix[i,i] == 1){
      index<-union(index,c(i))
    }else{
      index<-index}
  }
  reState<-(colnames(trMatrix)[index])

  if(is.null(index)){
    cat("No absorbing state exists in this Markov chain.")
    cat("\n")
  }else{cat("Absorbing state(s) : ")
    cat(reState,sep=", ")
    cat("\n")
  }
}

##plot method
getMC <- function(trMatrix, round = FALSE) {

  trMat <- trMatrix*500
  if(round == TRUE) {
    trMat <- round(trMat, 2)
  }

  MC <- igraph::graph.adjacency(adjmatrix = trMat, weighted = TRUE, mode = "directed")
  return(MC)
}

#'@import igraph
#'@export
plot.MarkovChain<-function(x,...,main,package="igraph"){

  trMat <- getMC(trMatrix = x, round = TRUE)
  Label <- round(igraph::E(trMat)$weight / 500, 2)

  igraph::plot.igraph(x = trMat, edge.label = Label,edge.arrow.size = 0.3,edge.arrow.width = 0.8,main = main)
}

##print method
#'
#'@export
print.MarkovChain<- function(x,...) {

  print(x[,])
  cat("\n")
}

## summary method
#'
#'@export
summary.MarkovChain<-function(object, ...){

  cat("States of this Markov chain : ")
  cat(displayState(object),sep=",")
  cat("\n")
  cat("\n")

  cat(transient(object))
  cat("\n")

  cat(recurrent(object))
  cat("\n")

  cat(absorbing(object))
  cat("\n")

  cat("Stationary distribution : ")
  cat("\n")
  return(stationary(object))
}
