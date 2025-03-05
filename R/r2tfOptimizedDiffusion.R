
##Set up parelell machine
##forkCluster <- makeCluster(32,'FORK',timeout = 100000)
##registerDoParallel(forkCluster)
##stopCluster(forkCluster)




#network <- data.frame(from = c('R','T1','T1','T3','T2','I1'),to = c('T1','T2','T3','T2','TF1','T3'))
#nodeW <- c(8,8,8,8,8,8)
#names(nodeW) <- c('R','T1','T2','T3','I1','TF1')
#signal <- c(1,0,0,0,0,0)

###Functions with feedback
#' Initialize a Signaling Network
#'
#' This function initializes a signaling network by assigning edge weights based on the product of connected node weights.
#'
#' @param network A dataframe with columns `from` and `to`, specifying network edges.
#' @param nodeW A named numeric vector of gene expression levels for each node.
#' 
#' @return A list containing:
#'   - `network`: A `data.table` object with edges defined by columns from and to.
#'   - `signal`: A named numeric vector typically initialized to zero. except for receptor
#'
#' @export
#' @import data.table
#' @examples
#' network <- data.frame(from = c('A', 'B'), to = c('B', 'C'))
#' nodeW <- c(A = 1, B = 2, C = 3)
#' initializeSignallingNetwork(network, nodeW)

initializeSignallingNetwork <- function(network,nodeW){

  if(!all(c(network$from,network$to) %in% names(nodeW))){
    stop('network genes missing in expression (nodeW)')
  }

  for(ii in 1:nrow(network)){
    network$edgeW[ii] <- nodeW[network$from[ii]]*nodeW[network$to[ii]]
  }

  ##Consider this one
  network$edgeW <- network$edgeW/(100*(max(network$edgeW)))
  network$flux <- 0


  signal <- rep(0,length(nodeW))
  names(signal) <- names(nodeW)
  nodes <- data.frame(nodeW= nodeW,signal = signal)
  network <- data.table(network)
  return(list(network= network,signal = signal))
}


#' Simulate Signal Flow in a Feedback Network
#'
#' This function simulates the propagation of a signal through a network with feedback regulation.
#'
#' @param network A list output from `initializeSignallingNetwork()`.
#' @param outputNode A character vector specifying the output nodes (TFs).
#' @param inputNode A character vector specifying the input nodes (e.g., receptors).
#' @param feedbackNodes A character vector specifying nodes that act as negative regulators.
#' @param inputSignal Numeric, amount of initial signal (default 0.99).
#' @param n Integer, number of iterations for signal diffusion (default 2000).
#'
#' @return A list containing:
#'   - `Enode`: A matrix tracking signal changes over time at output nodes.
#'   - `signal`: The final state of all nodes after diffusion.
#'
#' @import foreach
#' @import doParallel
#' @import data.table
#' @export
signalOnFeedbackNetwork <-function(network,outputNode,inputNode,
				   feedbackNodes = NULL ,
				   inputSignal = 0.99, n = 2000){



  signal <- network$signal
  network <- copy(network$network)
  signal[inputNode] <- inputSignal
  setkey(network,'from')
  Enode <- matrix(0,nrow = n,ncol = length(outputNode))
  colnames(Enode) <- outputNode

  for(ii in 1:n){
    network$delta <- signal[network$from]-signal[network$to]
    network$flux <- network$delta*network$edgeW
    f1 <- network[,.(flux = sum(flux)),by = from]
    r1 <- network[,.(flux = sum(flux)),by = to]
    signal[f1$from] <- signal[f1$from] - f1$flux
    signal[r1$to] <- signal[r1$to]+r1$flux
    if(!is.null(feedbackNodes)){
        signal[feedbackNodes] <- 0
    }
    Enode[ii,outputNode] <- signal[outputNode]
  }
  return(list(Enode=t(Enode),signal = signal))
}

#' Compute Diffusion Map for Signal Propagation
#'
#' Note on Parallel Processing in Network Diffusion
#'
#' This function runs a parallelized signal propagation model using `foreach` for 
#' efficient computation.
#'
#' @note This function requires a parallel backend to be registered for `foreach` to work.
#' Users can register a parallel backend using:
#' 
#' \code{
#' library(doParallel)
#' cl <- makeCluster(detectCores() - 1)  # Use available cores
#' registerDoParallel(cl)
#' }
#' The code is parallelized by sample so more cores than you have samples in the
#' will give no benefit it may therefore be best to dial down the
#' default on large clusters 
#'
#' @param receptors Character vector specifying receptor nodes.
#' @param TFs Character vector specifying transcription factor nodes.
#' @param feedbackNodes Character vector specifying nodes that act as signal sinks.
#' @param M A numeric matrix of gene expression data (rows = genes, cols = samples).
#' @param network A data frame representing the network topology with columns from and to.
#' @param inputSignal Numeric, amount of initial signal (default 0.99).
#' @param nTicks Integer, number of time steps for signal propagation (default 2000).
#' @param x Numeric vector, weights for network edges.
#'
#' @return A list containing:
#'   - `diffCurves`: Signal propagation curves for each sample.
#'   - `rowAnnotation`: Metadata linking curves to samples, TFs, and receptors.
#'
#' @import foreach
#' @import doParallel
#' @import data.table
#' @export
diffusionMap2 <- function(receptors, TFs,feedbackNodes, M, network,
			  inputSignal = 0.99, nTicks=2000,x = rep(1,nrow(network))){
  setDTthreads(1)
  if(!all(receptors %in% rownames(M))) {
    stop ("Please include receptors in the gene expression data (M) ")
  }
  if(!all(TFs %in% rownames(M))){
    stop("Please include TFs in the gene expression data (M) ")
  }
  if(!all(c(network[,1], network[,2]) %in% rownames(M))) {
    stop("Please include all network nodes in the gene expression data (M) ")
  }

 

  testf <- function(sampleName){
    # print(SignalOnNetwork)
    nodeW <- M[unique(c(network[,1],network[,2])),sampleName]
    startingNet <- initializeSignallingNetwork(network,nodeW)
    startingNet$network$edgeW <- startingNet$network$edgeW*x   
    aa <-   signalOnFeedbackNetwork(startingNet,outputNode = TFs,
                                inputNode = receptor,
                                feedbackNodes = feedbackNodes,
                                inputSignal = inputSignal,n = nTicks)
    diff_time <- aa[[1]]
  }
   sample_dtime <- vector('list',ncol(M))
   ## names(sample_dtime) <- colnames(M)
  
    ##Initialize storage array
    diffCurves <- NULL
    rowAnnotation <- data.frame(sampleid= NULL,TF = NULL, receptor = NULL)
    
    for(receptor in receptors ){
      diff_time <- foreach(sampleName = colnames(M), 
			   .export=c('initializeSignallingNetwork',
				     'signalOnFeedbackNetwork',
				     'setkey',
				     'data.table'),
			   .combine = rbind) %dopar% testf(sampleName)
     ## print((diff_time))
      diffCurves <- rbind(diffCurves,diff_time)
      tempAnnotation <- data.frame(sampleid = rep(colnames(M),
				   each = length(TFs)),
                                   TF = rownames(diff_time),
				   receptor = rep(receptor,
				   times = nrow(diff_time)))
      rowAnnotation <- rbind(rowAnnotation,tempAnnotation)
    }
    
    return(list(diffCurves = diffCurves,rowAnnotation = rowAnnotation))
   } 
    
#' Extract Diffusion Curves for a Specific Receptor-TF Pair
#'
#' This function retrieves diffusion curves for a specified receptor-transcription factor (TF) pair
#' from a previously computed diffusion analysis.
#'
#' @param X A list containing `diffCurves` (signal propagation data) and `rowAnnotation` (sample metadata).
#' @param receptor Character, the receptor of interest (default: `"LY96"`).
#' @param TF Character, the transcription factor of interest (default: `"NFKB1"`).
#' @param samples Character vector specifying sample IDs to extract (default: `NULL`, which includes all available samples).
#'
#' @return A list with:
#'   - `diffCurves`: A subset of the diffusion curves matrix corresponding to the specified receptor-TF pair and samples.
#'   - `rowAnnotation`: A subset of the row annotation dataframe for the selected samples.
#'
#' @export
#' @examples
#' # Example assuming `X` is a list with precomputed diffusion curves
#' selected_curves <- getCurves(X, receptor = "TLR4", TF = "RELA", samples = c("Sample1", "Sample2"))
#' str(selected_curves)  



  getCurves <- function(X,receptor= 'LY96',TF = 'NFKB1',samples = NULL){
      if(is.null(samples)){
          samples <- unique(X$rowAnnotation[,1])
      }
      
      ##Check if inputs exist and sample ids are unique
      
      ind <- which(X$rowAnnotation$sampleid %in% samples &
                    X$rowAnnotation$TF == TF &
                    X$rowAnnotation$receptor == receptor)
      
      Y <- list(diffCurves = X$diffCurves[ind,],
                rowAnnotation = X$rowAnnotation[ind,])
      rownames(Y$diffCurves) <- Y$rowAnnotation$sampleid          
      rownames(Y$rowAnnotation) <- Y$rowAnnotation$sampleid
  
      return(Y)
 }     
  
     ##Make sure data is returned in the order 
#' Summarize Diffusion Model Performance
#'
#' This function evaluates the relationship between network diffusion results and observed transcription factor (TF) activation data.
#' It computes the correlation, slope, and statistical significance (p-value) for each TF-receptor pair.
#'
#' @param res A list containing diffusion results, output from `diffusionMap2()`.
#' @param predData A numeric matrix with observed transcription factor activation levels (rows = samples, cols = TFs).
#' @param TFs A character vector of transcription factors to evaluate.
#' @param receptor A character string specifying the receptor of interest.
#'
#' @return A dataframe with rows corresponding to TFs and columns:
#'   - `rSq`: Pearson correlation coefficient between network diffusion results and observed TF activation.
#'   - `slope`: The slope of the linear regression fit.
#'   - `pValue`: The p-value from the linear regression analysis.
#'
#' @importFrom stats cor lm summary.lm var
#' @export
#' @examples
#' # Example assuming `res` contains diffusion results and `predData` is a TF activity matrix
#' TFs <- c("RELA", "NFKB1", "JUN")
#' receptor <- "TLR4"
#' summary_table <- summarizeDiffusionResult(res, predData, TFs, receptor)
#' print(summary_table)
summarizeDiffusionResult <- function(res,predData,TFs,receptor){

        ##Report correlation
        ##P value
        ##Report slope
        resTable <- data.frame(
                              rSq = rep(0,length(TFs)),
                               slope = rep(0,length(TFs)),
                               pValue = rep(0,length(TFs)),
                row.names = TFs)



        for(TF in TFs){
            XI <- getCurves(X = res,TF = TF,receptor = receptor)
            xiIQR <- apply(XI$diffCurves,2,var)
            ind <- which.max(xiIQR)
            fit <- lm(predData[,TF]~XI$diffCurves[,ind])
            fit <- summary(fit)

            resTable[TF,'rSq'] <- cor(XI$diffCurves[,ind],predData[,TF])
            if(!is.na(cor(XI$diffCurves[,ind],predData[,TF]))){
                resTable[TF,'slope'] <- fit$coefficients[2,1]
                resTable[TF,'pValue'] <- fit$coefficients[2,4]
            }
        }
        return(resTable)
}

#' Evaluate Network Connectivity in an Optimizer
#'
#' This function computes the performance score of a signaling network 
#' by comparing predicted transcription factor (TF) activation to observed data.
#' It is designed for use in optimization routines.
#'
#' @param x Numeric vector of edge weight modifiers for network optimization.Typical range 0 to 2.
#' @param network Data frame with two columns (`from`, `to`), representing network topology.
#' @param receptor Character string specifying the receptor gene ID.
#' @param TFs Character vector of transcription factor gene IDs.
#' @param predData Numeric matrix of observed TF activation (rows = samples, cols = TFs).
#' @param M Numeric matrix of gene expression data (rows = genes, cols = samples).
#' @param feedbackNodes Character vector of gene IDs acting as signal sinks (non-transmitting feedback nodes).
#' @param inputSignal Numeric, amount of initial signal at the receptor (default: `1`).
#' @param nTicks Integer, number of time steps for signal propagation (default: `2000`).
#' @param alpha Numeric, regularization weight for network sparsity.
#'
#' @return A numeric score representing the network's performance:
#'   - Lower values indicate better agreement with observed data.
#'   - Higher values indicate poorer agreement.
#'
#' @export
evaluateNetwork2 <- function(x,network,receptor,TFs,predData,M,
                            feedbackNodes,inputSignal,nTicks,alpha){
    ##  initialize network
    ##  x: network edge weight modifiers
    ##  network: network data frametwo columns from and to
    ##  receptor: receptor geen ID
    ##  TFs : TF gene ID
    ##  predData: vector number of samples
    ##  M: gene expression data
    ## feedbackNodes: vector of gene IDs signal sinks not transmitters used to feedback mechanisms
    ## inputSignal: amount of signal on node usually 1
    ##  
    ##nTicks : number of timeticks deafault to 2000 for now
    #print(rownames(M))

    ## Assume one or more TFs coming in

    res <- diffusionMap2(receptors = receptor,
              TFs = TFs,
              feedbackNodes= feedbackNodes,
              M = M,
              network = network,
              inputSignal = inputSignal,
              nTicks= nTicks,
              x = x)

    TFscore <- rep(0,length(TFs))
    names(TFscore) <- TFs

    for (TF in TFs){
        XI <- getCurves(X = res,TF = TF,receptor = receptor)
        xiIQR <- apply(XI[[1]],2,var)
        ind <- which.max(xiIQR)
        PC <- XI[[1]][,ind]
     ##Figure out multiple correlations
     ##Stick with one receptor
        if(var(PC) > 0){
           TFscore[TF] <- cor(PC,predData[,TF])
        }
        if(var(PC) == 0){
           TFscore[TF] <- -0.999
        }
        if(is.na(TFscore[TF])){
                TFscore[TF] <- -0.999
        }
        TFscore[TF] <- atanh(TFscore[TF])
     }

     score <- 1-tanh(mean(TFscore)) ## Should not use mean of correlations

     score <- mean(abs(x))*alpha + (1-alpha)*score
     return(score)
}

#' Optimize Network Edge Weights to Improve TF Activation Prediction
#'
#' This function optimizes the edge weights of a signaling network to maximize the correlation
#' between network connectivity and transcription factor (TF) activation. It first runs the 
#' raw diffusion model, then optimizes the network using the `nloptr` package, and finally 
#' evaluates the optimized network on a test dataset if provided.
#'
#' Note on Parallel Processing in Network Diffusion
#'
#' This function runs a parallelized signal propagation model using `foreach` for 
#' efficient computation.
#'
#' @note This function requires a parallel backend to be registered for `foreach` to work.
#' Users can register a parallel backend using:
#' 
#' \code{
#' library(doParallel)
#' cl <- makeCluster(detectCores() - 1)  # Use available cores
#' registerDoParallel(cl)
#' }
#' The code is parallelized by sample so more cores than you have samples in the
#' (ex vivo) training data will give no benefit it may therefore be best to dial down the
#' default on large clusters 
#' 
#' To stop parallel execution:
#' 
#' \code{
#' stopCluster(cl)
#' }
#'
#' @param network A dataframe with two columns (`from`, `to`), representing the signaling network.
#' @param receptor Character, the receptor gene ID.
#' @param TFs Character vector, transcription factor gene IDs.
#' @param feedbackNodes Character vector, gene IDs acting as feedback regulators (signal sinks).
#' @param Mtrain A numeric matrix of gene expression data (rows = genes, cols = training samples).
#' @param tfAct A numeric matrix of observed TF activation levels (rows = training samples, cols = TFs).
#' @param alpha Numeric, regularization weight for network sparsity (default: `0.01`).
#' @param Mtest A numeric matrix of gene expression data for the test dataset (optional, default: `NULL`).
#' @param patientData Optional patient dataset for validation (default: `NULL`).
#' @param maxeval Integer, maximum number of function evaluations in the optimization process (default: `250`).
#' @param nTicks Integer, number of time steps for signal propagation in diffusion model (default: `2000`).
#' @param x0 Numeric vector, initial values for network edge weight optimization (default: `rep(1, nrow(network))`).
#'
#' @return A list containing:
#'   - `opt`: Optimization results from `nloptr()`.
#'   - `corRaw`: Correlation results before optimization.
#'   - `corOpt`: Correlation results after optimization.
#'   - `rawRun`: Results from initial (unoptimized) network diffusion.
#'   - `optRun`: Results from optimized network diffusion.
#'   - `testSummary`: A list containing raw and optimized network connectivity dataset results for each TF (if `Mtest` is provided).
#'
#' @export
#'
#' @examples
#' # Example using example dataset (assuming r2tfData contains required inputs)
#' data(r2tfData)
#' optimized_results <- optimizeR2TF(
#'     network = r2tfData$network,
#'     receptor = r2tfData$receptors,
#'     TFs = r2tfData$TFs,
#'     feedbackNodes = r2tfData$feedbackNodes,
#'     Mtrain = r2tfData$M,
#'     tfAct = r2tfData$predData
#' )
#' print(optimized_results)
optimizeR2TF<- function(network,
                       receptor,
                       TFs,
                       feedbackNodes,
                       Mtrain,
                       tfAct,
                       alpha = 0.01,
                       Mtest = NULL,
                       patientData = NULL,
                       maxeval = 250,
                       nTicks = 2000,
		       x0 = NULL){
        #run raw diffusion
        if(is.null(x0)){
	    x0 <- rep(1,nrow(network))	   
        }
        opts <- list(algorithm = 'NLOPT_LN_BOBYQA' ,
                         maxeval = maxeval)

        rawRun <- diffusionMap2(receptors = receptor, TFs = TFs,feedbackNodes = feedbackNodes,
             M = Mtrain, network = network,inputSignal = 0.99, nTicks=nTicks)
        

        corRaw <- summarizeDiffusionResult(res =rawRun,
                         predData = tfAct,
                         TFs = TFs,
                         receptor = receptor)

        optRes <- nloptr::nloptr(x0 = x0,eval_f = evaluateNetwork2,
                             ub = rep(2,nrow(network)), 
			     lb = rep(0,nrow(network)),
                             opts = opts,
                             network = network,
                             receptor = receptor,
                             TFs = TFs,
                             predData = tfAct,
                             M = Mtrain,
                             feedbackNodes = feedbackNodes,
                             inputSignal = 0.99,
                             nTicks = nTicks,
                             alpha = alpha)

         optRun <- diffusionMap2(receptors = receptor,
                             TFs = TFs,
                             feedbackNodes= feedbackNodes,
                             M = Mtrain,
                             network = network,
                             inputSignal = 0.99,
                             nTicks= nTicks,
                             x = optRes$solution)

         corOpt <- summarizeDiffusionResult(res =optRun,
                         predData = tfAct,
                         TFs =TFs,
                         receptor = receptor)

                testSummary <- vector('list',length(TFs))
          names(testSummary) <- TFs

    if(!is.null(Mtest)){

              Hraw <- diffusionMap2(receptors = receptor , TFs = TFs,
                      feedbackNodes = feedbackNodes,
                      M = Mtest, network = network,
                      inputSignal = 0.99, nTicks=nTicks)


              Hopt <- diffusionMap2(receptors = receptor , TFs = TFs,
                      feedbackNodes = feedbackNodes,
                      M = Mtest, network = network,
                      inputSignal = 0.99, nTicks=nTicks,
                      x = optRes$solution)

              for (TF in TFs){
                  Xraw <- getCurves(X = Hraw,TF = TF,receptor = receptor)
                  xVar <- apply(Xraw[[1]],2,var )
                  tmax <- which.max(xVar)
                  Sraw <- Xraw[[1]][,tmax]

                  Xopt <- getCurves(X = Hopt,TF = TF,receptor = receptor)
                  xVar <- apply(Xopt[[1]],2,var )
                  tmax <- which.max(xVar)
                  Sopt <- Xopt[[1]][,tmax]
    
                  testSummary[[TF]] <- list(Sraw = Sraw,Sopt = Sopt)
              }

        } ##If test
  
   return(list(opt = optRes,corRaw = corRaw, corOpt = corOpt,
	       rawRun = rawRun,optRun = optRun, 
	       testSummary = testSummary))
}



