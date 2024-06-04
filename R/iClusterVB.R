mydir <- "./src/"
sourceCpp(paste0(mydir,"CAVI_algorithm_standard_version.cpp"))
sourceCpp(paste0(mydir,"CAVI_algorithm_global_version.cpp"))

CAVI_algorithm_global <- cmpfun(CAVI_algorithm_global)
CAVI_algorithm_standard <- cmpfun(CAVI_algorithm_standard)

#---------------------------------------------------------------------------#
iClusterVB <- function(
    mydata, 			 # input data - list of length R data set (each data set is of N times p_r dimensional)
    dist,				 # type of data, or distribution: dist = "continuous", "multinomial", "poisson" (vector of length R)
    K = 10, 			 # maximum number of clusters
    initial_method = "VarSelLCM",    # initialization method: VarSelLCM, random, k-prototype, k-means, mclust, lca,
    VS_method = 0, 			 # variable selection method: 0 = clustering no variable selection, 1 = clustering with variable selection
    initial_cluster = NULL,          # initial cluster membership, if it is not NULL, it will overwrite the previous initial values setting for this parameter
    initial_vs_prob = NULL,          # initial variable selection probability, a scalar
    initial_fit = NULL,              # initial values based on previously fitted iClusterVB model (an iClusterVB object)
    initial_omega = NULL, 		 # constomized initial values for variable inclusion probabilities, if it is not NULL, it will overwrite the previous initial values setting for this parameter
    # if VS_method  = 1, initial_omega is a list with length R, each element of the list is an array with dim=c(N,p[[r]])), r = 1,...,R
    input_hyper_parameters = NULL,	 # input hyper-parameters of the prior distributions for the model
    max_iter = 200,			 # maximum iteration of the VB algorithm
    early_stop = 1,			 # whether stop the algorithm when converge or continue until reach max_iter: 1 - algorithm stops when converges,
    # 0 - algorithm stops when reaches the maximum iteration (regardless of converges or not)
    per = 10,                        # print information every per iteration
    convergence_threshold = 1e-4	   # Define a convergence threshold for the change in ELBO
){

  test <- 0
  if (test==1){
    K = 10
    initial_method = "VarSelLCM"
    VS_method = 0
    initial_cluster = NULL
    initial_vs_prob = NULL
    initial_fit = NULL
    initial_omega = NULL
    input_hyper_parameters = NULL
    max_iter = 200
    early_stop = 1
    per = 10
    convergence_threshold = 1e-4
  }

  cat(paste(rep('-',60),sep='',collapse=''), '\n');
  cat("Pre-processing and initilizing the model",  '\n')
  cat(paste(rep('-',60),sep='',collapse=''), '\n');

  # evaluating if there is missing data
  if (sum(is.na(do.call(cbind, mydata))) > 0 ) {
    stop("Error: missing data is not allowed\n")
  }

  if (initial_method == "VarSelLCM" | initial_method == "kproto"){mydata_new <- mydata}

  #------
  R <- length(mydata)	# number of dataset
  # function evaluating if the number of categories for the categorical variables are identical or not
  is_all_equal <- function(lst) {
    all(sapply(lst, function(x) identical(unlist(lst[[1]]), unlist(x))))
  }
  # vector store the number of categories for each categorial variable (for each dataset)
  num.categories <- list()
  for (r in 1:R){
    if(dist[r] == "multinomial"){
      indx <- which(which(dist == "multinomial")==r)
      num.categories[[indx]] <- apply(mydata[[r]], 2, function(x) length(unique(x)));
      if (initial_method == "VarSelLCM" | initial_method == "kproto") mydata_new[[r]] <- as.data.frame(lapply(data.frame(mydata_new[[r]]), as.factor))   # convert to categorical data to factor for further analysis
    }
  }
  M <- lapply(num.categories, function(x) x[1]); M # M is a list

  # Evaluating if the sample size of the input datasets are the same
  sample.size <- lapply(mydata, function(x) dim(x)[1])   # sample size each dataset (in a list)
  if (is_all_equal(sample.size) == "FALSE") {
    stop("Error: sample size of the input datasets are not the same \n")
  }

  if (initial_method == "VarSelLCM" | initial_method == "kproto"  ){
    mydata_combine <- do.call(cbind, mydata_new)   # categorical data have been converted to factor
    mydata_combine <- as.data.frame(mydata_combine)
    colnames(mydata_combine) <- paste0("x",1:dim(mydata_combine)[2])} else{mydata_combine <- do.call(cbind, mydata)}

  N <- sample.size[[1]]; N

  # dimension of each dataset
  p <- lapply(mydata, function(x) dim(x)[2]);	p		# dimension of each dataset (in a list)

  #---------------------------------------------------------------------------#
  p_total <- sum(unlist(p)); p_total					# total dimension of all datasets
  p_continuous_total <- sum(dist=="gaussian")
  p_categorical_total <- sum(dist=="multinomial")
  p_count_total <- sum(dist=="poisson")


  if ( (initial_method == "lca") & (p_continuous_total + p_count_total > 0 ) ) {
    stop("lca method can only be used as an initialization method when there are only categorical datasets \n")
  }

  data_summary <- list(
    N = N,
    M = M,
    p = p,
    R = R,
    p_total = p_total,
    p_continuous_total = p_continuous_total,
    p_categorical_total = p_categorical_total,
    p_count_total = p_count_total
  )

  #---------------------------------------------------------------------------#
  # Hyper-parameters
  #---------------------------------------------------------------------------#
  if (is.null(input_hyper_parameters) == TRUE ){
    # default hyper-parameter setting
    input_hyper_parameters <- list(
      alpha_00 = 0.001,		# hyper-parameter for the mixing proportion
      mu_00 = 0,			# mean hyper-parameter for the mean of the Gaussian distribution
      s2_00 = 100,			# variance hyper-parameter for the mean of  the Gaussian distribution
      a_00 = 1,			# shape parameter for Iv-Gamma for the variance of the Gaussian distribution
      b_00 = 1,			# scale parameter for Iv-Gamma for the variance of the Gaussian distribution
      kappa_00 = 1, 		# hyper-parameter for the Dirichlet distribution, for the categorical outcomes
      u_00 = 1,			# shape parameter for Iv-Gamma for the variance of the Gaussian distribution
      v_00 = 1)			# scale parameter for Iv-Gamma for the rate parameter of the Poisson distribution
  }

  # Hyper-parameters for dataset with continuous variables
  mu_0 <- list()  # Prior mean for cluster-specific means
  s2_0 <-  numeric()		# Prior variance for cluster-specific means
  a_0 <-  numeric() 		# Shape parameter for inverse-gamma prior for variances
  b_0 <-  numeric() 		# Rate parameter for inverse-gamma prior for variances

  # Hyper-parameters for dataset with categorical variables (Dirichlet prior distribution)
  kappa_0 <- list()

  # Hyper-parameters for dataset with count variables (Gamma prior distribution)
  u_0 <- numeric()			# hyper-parameter of the gamma distribution for the rate parameter
  v_0 <- numeric()		    	# hyper-parameter of the gamma distribution for the rate parameter

  alpha_0 <- rep(input_hyper_parameters$alpha_00,K)

  for (r in 1:R){
    #------------------------#
    if(dist[r] == "gaussian"){

      indx <- which(which(dist == "gaussian")==r)	# the sequential order of the guassian distributed dataset
      # for dataset with continuous variables: specify prior for means and variances
      mu_0[[indx]] <- rep(input_hyper_parameters$mu_00, p[[r]])
      s2_0[indx]  <- input_hyper_parameters$s2_00
      a_0[indx] <- input_hyper_parameters$a_00
      b_0[indx] <- input_hyper_parameters$b_00
    }
    # for dataset with categorical variables: specify prior (Dirichlet distribution) for the probability of categories of each variable
    if(dist[r] == "multinomial"){
      indx <- which(which(dist == "multinomial")==r)
      kappa_0[[indx]] <- array(input_hyper_parameters$kappa_00, dim = c(K, p[[r]], M[[indx]])) # Hyperparameters (Dirichlet distribution) for categories
    }
    if(dist[r] == "poisson"){
      # for dataset with count variables: specify prior (gamma prior) for the rate parameter of each variable
      indx <- which(which(dist == "poisson")==r)
      u_0[indx] <- input_hyper_parameters$u_00
      v_0[indx] <- input_hyper_parameters$v_00
    }
  }

  #---------------------------------------------------------------------------#
  hyper_parameters <- list(
    mu_0 = mu_0,
    s2_0 = s2_0,
    a_0 = a_0,
    b_0 = b_0,
    kappa_0 = kappa_0,
    u_0 = u_0,
    v_0 = v_0,
    alpha_0 = alpha_0 )

  #---------------------------------------------------------------------------#
  # Define a small value to add for numerical stability
  epsilon <- 1e-10

  #---------------------------------------------------------------------------#
  # Initialize cluster allocation
  #---------------------------------------------------------------------------#
  if (is.null(initial_cluster) == TRUE & is.null(initial_fit)==TRUE) {

    #---------------------------------------------------------------------------#
    # Initialize cluster allocation (random initialization)
    #---------------------------------------------------------------------------#
    if (initial_method == "random"){
      # initial values related to the clustering
      zz <- initial_cluster <- sample(1:K, N, replace=TRUE)			# random sample of the cluster allocation
    }
    #---------------------------------------------------------------------------#
    # Initialize cluster allocation (VarSelLCM initialization)
    #---------------------------------------------------------------------------#
    if (initial_method == "VarSelLCM" & VS_method == 0){
      #fit_VarSelLCM <- suppressWarnings(VarSelCluster(mydata_combine, gvals = 1:K, vbleSelec = FALSE, crit.varsel = "BIC"))
      fit_VarSelLCM <- suppressWarnings(VarSelCluster(mydata_combine, gvals = K, vbleSelec = FALSE, crit.varsel = "BIC"))
      zz <- initial_cluster <- as.numeric(fitted(fit_VarSelLCM));  table(zz)
    }
    if (initial_method == "VarSelLCM" & VS_method != 0){
      #fit_VarSelLCM <- suppressWarnings(VarSelCluster(mydata_combine, gvals = 1:K, vbleSelec = TRUE, crit.varsel = "BIC", iterKeep = 20))
      fit_VarSelLCM <- suppressWarnings(VarSelCluster(mydata_combine, gvals = K, vbleSelec = FALSE, crit.varsel = "BIC", iterKeep = 50))
      zz <- initial_cluster <- as.numeric(fitted(fit_VarSelLCM));  table(zz)
      #omega_VarSelLCM <-  as.numeric(fit_VarSelLCM@model@omega)
    }
    #----------------------------------------------------------------------------#
    # Initialize cluster allocation (K-Medoids)
    #----------------------------------------------------------------------------#
    if (initial_method == "PAM"){
      gower_dist <- daisy(mydata_combine,metric = "gower")
      fit.PAM <- pam(gower_dist, diss = TRUE, k = K, cluster.only=TRUE)
      zz <- initial_cluster <- as.numeric(fit.PAM$clustering); table(zz)
    }
    #----------------------------------------------------------------------------#
    # Initialize cluster allocation (kproto initialization)
    #----------------------------------------------------------------------------#
    if (initial_method == "kproto"){
      # install.packages("clustMixType")
      fit.kproto <- kproto_gower(x = data.frame(mydata_combine), k  = K, lambda = rep(1,p_total));
      zz <- initial_cluster <- as.numeric(fit.kproto$cluster); table(fit.kproto$cluster)
    }
    #----------------------------------------------------------------------------#
    # Initialize cluster allocation: continuous data only (k-means initialization)
    #----------------------------------------------------------------------------#
    #if (initial_method == "kmeans" & (p_categorical_total) == 0){
    if (initial_method == "kmeans"){
      fit.km <- kmeans(data.frame(mydata_combine), centers = K)
      zz <- initial_cluster <- as.numeric(fit.km$cluster)
    }
    #----------------------------------------------------------------------------#
    # Initialize cluster allocation: continuous data only (mclust initialization)
    #----------------------------------------------------------------------------#
    #if (initial_method == "mclust" & (p_categorical_total) == 0){
    if (initial_method == "mclust"){
      fit.mclust <- Mclust(as.matrix(mydata_combine),G = K, verbose = FALSE) ;
      zz <- initial_cluster <- as.numeric(fit.mclust$classification)
    }
    #------------------------------------------------------------------------------------#
    # Initialize cluster allocation: categorical data only	(lca initialization)
    #------------------------------------------------------------------------------------#
    if (initial_method == "lca" & (p_continuous_total + p_count_total) == 0){
      # recode 0 to a category with positive number
      mydata_combine[mydata_combine==0] <- max(mydata_combine) + 1
      fit.lca <- poLCA(as.matrix(mydata_combine) ~ 1, data= NULL, nclass= K , verbose = FALSE, calc.se = FALSE)
      zz <- initial_cluster <- as.numeric(fit.lca$predclass); table(zz)
    }
  }
  if (is.null(initial_cluster) == FALSE) { zz <- as.numeric(initial_cluster)}
  if (is.null(initial_fit) == FALSE) {zz <- as.numeric(initial_fit$cluster)} # over-write the previous class assignment

  #------------------------------------------------------------------------------------#
  # Initialize other cluster allocation parameters
  #------------------------------------------------------------------------------------#
  zz <- factor(zz, levels=1:K)
  phi <- as.matrix(model.matrix(~ zz - 1)) # Cluster assignment probabilities
  ppi <- apply(phi, 2 , mean)				    # Mixing proportions
  alpha <-  as.numeric(table(zz))
  log_phi <- log(phi + epsilon)
  #------------------------------------------------------------------------------------#
  # Initial values related to the variable selection
  #------------------------------------------------------------------------------------#
  omega <- list()			# similar to pi
  ss <- list()				# similar to zz
  rho <- list()
  omega_tmp <- list()
  e <- list()
  gamma <- list()
  xi <- list()
  #------
  if (VS_method == 1 & is.null(initial_fit) == TRUE){
    for (r in 1:R){
      if(is.null(initial_vs_prob) == TRUE) {
        #if(N >= p_total)omega[[r]] <- array(0.9, dim=c(N,p[[r]]))
        #if(N < p_total) omega[[r]] <- array(0.1, dim=c(N,p[[r]]))
        omega[[r]] <- array(0.5, dim=c(N,p[[r]]))
      }  # similar to phi
      if(is.null(initial_vs_prob) == FALSE){omega[[r]] <- array(initial_vs_prob, dim=c(N,p[[r]])) }  # similar to phi                                  }
      ss[[r]] <- array(0, dim=c(N,p[[r]]))  	 # similar to zz
      ss[[r]][omega[[r]] >=0.5] <- 1
      rho[[r]] <- t(apply(omega[[r]], MARGIN = 2, FUN = mean));rho
    }
    omega_tmp <- omega 	# initial value only
    e <- gamma <- rho		# initial value only (for variational parameters)
    xi <- omega			    # initial value only
  }

  #----------------------------------------------------------------------------#
  #----------------------------------------------------------------------------#
  mu <- list() 	 		          # Cluster-specific means
  sigma2 <- list()  			    # Cluster-specific variances
  mu_tilde <- list() 		      # Posterior mean update for mean parameter
  s2_tilde <- list() 		      # Posterior variance update for mean parameter
  a_tilde <- b_tilde <- list()	# parameter for inverse gamma distribution
  theta <-	list()			        # parameter for dirichlet distribution of categorical data
  kappa <- 	list()
  theta_e <-	list()
  lambda  <-	list()
  lambda_e  <-	list()
  u_tilde  <-	list()			    # parameter for poisson distribution
  v_tilde  <-	list()			    # parameter for poisson distribution
  #----------------------------------------------------------------------------#
  for (r in 1:R){
    # for continuous variables
    if(dist[r] == "gaussian"){
      indx <- which(which(dist == "gaussian")==r)	# the sequential order of the Gaussian distributed dataset
      mu[[indx]] <-  matrix(0, ncol=p[[r]], nrow=K)
      mu[[indx]][1:length(unique(zz)),] <-  as.matrix(aggregate(mydata[[r]],list(zz),mean)[,-1])               # Cluster-specific means
      sigma2[[indx]] <-  matrix(0, ncol=p[[r]], nrow=K)
      sigma2[[indx]][1:length(unique(zz)),] <- as.matrix(aggregate(mydata[[r]],list(zz),var)[,-1])              # Cluster-specific variances
      a_tilde[[indx]] <-   matrix(0.01, ncol=p[[r]], nrow=K)  		# initial value only
      b_tilde[[indx]] <-   matrix(0.01, ncol=p[[r]], nrow=K)		# initial value only
      mu_tilde[[indx]] <- mu[[indx]] 					    # initial value only
      s2_tilde[[indx]] <- sigma2[[indx]]  				# initial value only
    }
    # for categorical data
    if(dist[r] == "multinomial"){
      indx <- which(which(dist == "multinomial")==r)
      tmp <- array(0, dim=c(K,p[[r]],M[[indx]]))

      theta[[indx]] <-  array(0, dim = c(K, p[[r]], M[[indx]]))
      kappa[[indx]] <-  array(0, dim = c(K, p[[r]], M[[indx]]))
      theta_e[[indx]] <-  array(0, dim = c(K, p[[r]], M[[indx]]))

      for (k in 1:K) {
        for (j in 1:p[[r]]) {
          for (m in 1:M[[indx]]){
            kappa[[indx]][k, j, m] <- kappa_0[[indx]][k, j, m] + sum((mydata[[r]][,j]==m) * phi[,k])
            theta_e[[indx]][k, j, m] <-  sum(mydata[[r]][,j]==m)/N
          }
          theta[[indx]][k,j,] <- 	kappa[[indx]][k, j, ]/sum(kappa[[indx]][k,j,])
        }
      }
    }
    # for count data
    if(dist[r] == "poisson"){
      indx <- which(which(dist == "poisson")==r)
      lambda[[indx]] <-  matrix(0, ncol=p[[r]], nrow=K)
      lambda[[indx]][1:length(unique(zz)),] <- as.matrix(aggregate(mydata[[r]], list(zz), mean)[,-1]); lambda
      lambda_e[[indx]] <- lambda[[indx]] 	# initial value only
      u_tilde[[indx]] <- lambda[[indx]] 					# initial value only
      v_tilde[[indx]] <- lambda[[indx]]/lambda[[indx]]  		# initial value only
    }
  }

  #------------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------------#
  # obtain initial values from the previous iClusterVB model
  if (is.null(initial_fit) == FALSE ){

    if (!(class(initial_fit) %in% c("CAVI_algorithm_standard","CAVI_algorithm_global")) == TRUE){
      stop("initial_fit should be either 'CAVI_algorithm_standard' or 'CAVI_algorithm_global' class \n")}

    # clustering parameters
    zz <- initial_cluster <- initial_fit$cluster
    phi <- initial_fit$model_parameters$phi
    ppi <- initial_fit$model_parameters$ppi
    alpha <-  as.numeric(table(zz))
    log_phi <- log(phi + epsilon)

    # initial values for cluster-specific parameters (continuous variables)
    mu = initial_fit$model_parameters$mu
    sigma2 = initial_fit$model_parameters$sigma2
    mu_tilde = initial_fit$model_parameters$mu_tilde
    s2_tilde = initial_fit$model_parameters$s2_tilde
    a_tilde = initial_fit$model_parameters$a_tilde
    b_tilde = initial_fit$model_parameters$b_tilde
    # initial values for cluster-specific parameters (categorical variables)
    theta = initial_fit$model_parameters$theta
    theta_e = initial_fit$model_parameters$theta_e
    kappa = initial_fit$model_parameters$kappa
    # initial values for cluster-specific parameters (count variables)
    lambda = initial_fit$model_parameters$lambda
    lambda_e = initial_fit$model_parameters$lambda_e
    u_tilde = initial_fit$model_parameters$u_tilde
    v_tilde = initial_fit$model_parameters$v_tilde

    # variable selection parameters
    # -------
    if(VS_method == 1 & class(initial_fit) == "CAVI_algorithm_standard"){
      for (r in 1:R){
        if(is.null(initial_vs_prob) == TRUE){	omega[[r]] <- array(0.5, dim=c(N,p[[r]])) }  # similar to phi
        if(is.null(initial_vs_prob) == FALSE){omega[[r]] <- array(initial_vs_prob, dim=c(N,p[[r]])) }  # similar to phi
        ss[[r]] <- array(0, dim=c(N,p[[r]]))  	 # similar to zz
        ss[[r]][omega[[r]] >=0.5] <- 1
        rho[[r]] <- t(apply(omega[[r]], MARGIN = 2, FUN = mean));rho
      }
      omega_tmp <- omega
      e <- gamma <- rho
      xi <- omega
    }
    # -------
    if(VS_method ==1 & class(initial_fit) == "CAVI_algorithm_global"){
      omega <- initial_fit$model_parameters$omega
      ss <- initial_fit$model_parameters$ss
      rho <- initial_fit$model_parameters$rho
      #omega_tmp <- initial_fit$model_parameters$omega_tmp
      #xi <- initial_fit$model_parameters$xi
      #e <- initial_fit$model_parameters$e
      #gamma <- initial_fit$model_parameters$gamma
      omega_tmp <- omega
      e <- gamma <- rho
      xi <- omega
    }
  }

  # Customized initial values for cluster membership
  # If it is not NULL, it will overwrite the previous initial values setting for this parameter
  # if (is.null(initial_cluster) == FALSE ){ zz <- initial_cluster}

  # Customized initial values for variable inclusion probabilities
  # If it is not NULL, it will overwrite the previous initial values setting for this parameter

  if (is.null(initial_omega) == FALSE ){
    omega <- initial_omega
    #----------
    if(VS_method == 1){
      for (r in 1:R){
        ss[[r]] <- array(0, dim=c(N,p[[r]]))  	 # similar to zz
        ss[[r]][omega[[r]] >=0.5] <- 1
        rho[[r]] <- t(apply(omega[[r]], MARGIN = 2, FUN = mean));rho
      }
      omega_tmp <- omega
      e <- gamma <- rho
      xi <- omega
    }
  }

  #------------------------------------------------------------------------------------#
  initial_values <- list(
    initial_method = initial_method,
    # initial values related to the clustering
    phi = phi,
    ppi = ppi,
    alpha = alpha,
    zz = zz,
    # initial values related to the variable selection
    omega =  omega,
    ss =  ss,
    rho =  rho,
    omega_tmp = omega_tmp,
    xi = xi,
    e = e,
    gamma = gamma,
    # initial values for cluster-specific parameters (continuous variables)
    mu = mu,
    sigma2 = sigma2,
    mu_tilde = mu_tilde,
    s2_tilde = s2_tilde,
    a_tilde = a_tilde,
    b_tilde = b_tilde,
    # initial values for cluster-specific parameters (categorical variables)
    theta = theta,
    theta_e = theta_e,
    kappa = kappa,
    # initial values for cluster-specific parameters (count variables)
    lambda = lambda,
    lambda_e = lambda_e,
    u_tilde = u_tilde,
    v_tilde = v_tilde
  )

  # make sure each dataset is a matrix before pass to Rcpp function
  #mydata <- lapply(mydata, function(y) if(is.factor(y)) as.matrix(as.numeric(as.character(y))) else y)
  #str(mydata)

  cat(paste(rep('-',60),sep='',collapse=''), '\n');
  cat("Running the CAVI algorithm",  '\n')
  cat(paste(rep('-',60),sep='',collapse=''), '\n');

  #------------------------------------------------------------------------------------#
  # CAVI algorithm parameters
  #------------------------------------------------------------------------------------#
  elbo <- numeric(max_iter)

  #------------------------------------------------------------------------------------#
  # Run CAVI algorithm
  #------------------------------------------------------------------------------------#
  start_time <- Sys.time()  # Record start time
  if (VS_method == 0){
    res_CAVI <- CAVI_algorithm_standard(
      mydata,
      data_summary,
      hyper_parameters,
      initial_values,
      dist,
      K,
      max_iter,
      early_stop,
      per,
      epsilon,
      convergence_threshold)
    class(res_CAVI) <- "CAVI_algorithm_standard"

  }


  if (VS_method == 1){
    res_CAVI <- CAVI_algorithm_global(
      mydata,
      data_summary,
      hyper_parameters,
      initial_values,
      dist,
      K,
      max_iter,
      early_stop,
      per,
      epsilon,
      convergence_threshold)
    class(res_CAVI) <- "CAVI_algorithm_global"
  }


  # return results
  res_CAVI

}
