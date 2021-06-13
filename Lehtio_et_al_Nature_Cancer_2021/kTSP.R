library(switchBox)

# data needs to be imported
pairs.mat <- read.table("./data/Most_13_Pairs_ktsp.txt", header = T)
all.prots <- read.table("./data/All_Proteins_in_k13.txt")

rem.na <- function(vector.x){
    to.remove <- which(vector.x == "Tie-Vote")
    if(length(to.remove) > 0){
        vector.x <- vector.x[,-to.remove]
    }else{
        vector.x <- vector.x
    }
    numeric.vector <- as.numeric(apply(vector.x, 1, function(x) as.numeric(as.character(x))))
    names(numeric.vector) <- names(vector.x)
    most.voted <- names(sort(table(numeric.vector), decreasing = T)[1:3])
    if(unname(sort(table(numeric.vector), decreasing = T)[1] == unname(sort(table(numeric.vector), decreasing = T)[3]))){
        return("Uncls")
    }else if(unname(sort(table(numeric.vector), decreasing = T)[1] == unname(sort(table(numeric.vector), decreasing = T)[2]))){
        c1 <- sort(as.numeric(most.voted)[1:2], decreasing = F)[1]
        c2 <- sort(as.numeric(most.voted)[1:2], decreasing = F)[2]
        comp <- paste(c1, c2, sep = "-")
        if(comp %in% colnames(vector.x)){
            return(as.numeric(as.character(vector.x[,comp])))
        }else{
            return("Uncls")
        }
    }else{
        return(most.voted[1])
    }
}

ktsp.13 <- function(data){

    #check feature pairs overlap
    # get the samples higher than 50% coverage

    ss.classification <- lapply(1:ncol(data), function(x){
        one.sample <- data[, x, drop = FALSE]
        one.sample <- na.omit(one.sample)
        im.ratio <- unname(unlist(sapply(1:nrow(pairs.mat), function(p){
            pair1 <- strsplit(pairs.mat[p,1], split = "-")[[1]][1]
            pair2 <- strsplit(pairs.mat[p,1], split = "-")[[1]][2]
            ifelse(length(intersect(rownames(one.sample), c(pair1, pair2))) > 0, 1, 0)
        })))
        no.ratio <- sum(im.ratio) / nrow(pairs.mat)

        temp.pair <- pairs.mat
        temp.pair$Imp.Ratio <- im.ratio
        temp.pair <- temp.pair[temp.pair$Imp.Ratio == 1, ]
        temp.pair$p1 <- unname(unlist(sapply(1:nrow(temp.pair), function(p){strsplit(temp.pair[p,1], split = "-")[[1]][1]})))
        temp.pair$p2 <- unname(unlist(sapply(1:nrow(temp.pair), function(p){strsplit(temp.pair[p,1], split = "-")[[1]][2]})))

        if(no.ratio >= 0.5){
            #Get the overlapped pair features and prepare the classifier
            diff.prots <- setdiff(as.character(all.prots$Prot), unique(c(temp.pair$p1, temp.pair$p2)))

            rules <- data.frame(Pairs = as.numeric(),
                                pair = as.numeric(),
                                freq.Freq = as.numeric(),
                                p1 = as.numeric(),
                                p2 = as.numeric())

            classifiers <- list()
            k <- 0
            for(m in seq_len(5)){
                t = m + 1
                for(n in t:6){
                    k <- k +1
                    temp.pairs <- pairs.mat[pairs.mat$pair == paste(m,n,sep = "-"),]

                    temp.pairs$No.Id <- unlist(sapply(1:nrow(temp.pairs),function(x){
                        ifelse(sum(unlist(strsplit(temp.pairs[x,1], split = "-")) %in% diff.prots == TRUE) > 0, TRUE, FALSE)
                    }))

                    temp.pairs <- temp.pairs[temp.pairs$No.Id == "FALSE", ]

                    #if(nrow(temp.pairs) == 0 ){
                    #    stop("! There is no rule found for some clusters")
                    #    return(NA)
                    #}

                    temp.pairs$p1 <- unlist(sapply(1:nrow(temp.pairs), function(x){
                        unlist(strsplit(temp.pairs[x,1], split = "-"))[1]}))

                    temp.pairs$p2 <- unlist(sapply(1:nrow(temp.pairs), function(x){
                        unlist(strsplit(temp.pairs[x,1], split = "-"))[2]}))

                    rules <- rbind(rules, temp.pairs[,c(1,2,4,6,7)])

                    TSPs <- matrix(c(as.character(temp.pairs$p1), temp.pairs$p2),
                                   nrow=nrow(temp.pairs),
                                   dimnames=list(c(temp.pairs$Pairs),
                                                 c("gene1", "gene2")))

                    name <- paste(nrow(temp.pairs), "TSPS", sep = "")
                    score <- rep(1.01, nrow(temp.pairs))
                    tieVote <- factor(rep("both",nrow(temp.pairs)), levels = c("both", m, n))
                    names(tieVote) <- dimnames(TSPs)[[1]]

                    labels <- c(m,n)
                    #print(paste(paste(m, n , sep = "-"), name , sep = ":"))

                    cls <- list(name = name,
                                TSPs = TSPs,
                                score = score,
                                tieVote = tieVote,
                                labels = labels)

                    classifiers[[k]] <- cls
                }
            }
            res <- list(rules, classifiers)

            #results <- tryCatch(ktsp.cls.prep(x = , all.prots = all.prots),
            #                    error=function(err) NA)
            rules <- res[[1]]
            classifiers <- res[[2]]

            # run k-TSP classification
            # imput one sample
            proteins  <- as.character(unique(c(rules$p1, rules$p2)))

            one.sample <- data[proteins, x, drop = FALSE]

            not.in.data <- setdiff(proteins, rownames(one.sample))

            not.in.data.df <- data.frame(rep(NA, length(not.in.data)), row.names = not.in.data)
            colnames(not.in.data.df) <- colnames(one.sample)

            one.sample <- rbind(one.sample, not.in.data.df)
            one.sample <- replace(one.sample, is.na(one.sample), 1)

            classification <- lapply(1:15, function(y){
                #Counting the signed TSP votes
                classifier <- classifiers[[y]]
                ktspStatDefault <- SWAP.KTSP.Statistics(inputMat = as.matrix(one.sample),
                                                        classifier = classifier,
                                                        CombineFunc=sum)

                #Apply the classifier to the new set
                testPrediction <- SWAP.KTSP.Classify(as.matrix(one.sample), classifier)

                if(nrow(classifier$TSPs) - unname(ktspStatDefault$statistics) == unname(ktspStatDefault$statistics)){
                    testPrediction <- "Tie-Vote"
                    names(testPrediction) <- colnames(one.sample)
                }

                cls.df <- data.frame(testPrediction,
                                     False_Vote = nrow(classifier$TSPs) - unname(ktspStatDefault$statistics),
                                     True_Vote = unname(ktspStatDefault$statistics),
                                     Sample = names(testPrediction))

            })
            all.cls.df <- do.call("cbind", classification)
            cls.df <- all.cls.df[, seq(from =1, to = 60, by = 4)]
            colnames(cls.df) <- c("1-2", "1-3", "1-4", "1-5", "1-6",
                                  "2-3", "2-4", "2-5", "2-6", "3-4",
                                  "3-5", "3-6", "4-5", "4-6", "5-6")
            cls.df$FinalPred <- rem.na(vector.x = cls.df)

            f.vote.df <- all.cls.df[, seq(from =2, to = 60, by = 4)]
            colnames(f.vote.df) <- paste(colnames(cls.df)[1: ncol(cls.df) -1 ], "False Vote", sep = "-")

            t.vote.df <- all.cls.df[, seq(from =3, to = 60, by = 4)]
            colnames(t.vote.df) <- paste(colnames(cls.df)[1: ncol(cls.df) -1 ], "True Vote", sep = "-")

            all.df <- cbind(cls.df, t.vote.df, f.vote.df)
            all.df <- all.df[, c(1:16, 17, 32, 18, 33, 19, 34, 20, 35, 21, 36, 22, 37, 23, 38, 24, 39, 25, 40, 26, 41, 27, 42, 28, 43, 29, 44, 30, 45, 31, 46)]
            all.df$TotalPAirs <- sum(im.ratio)
            all.df$PairsCoverage <- no.ratio
            all.df

        }else{
            return(NA)
        }
    })
    cat("Classification done!\n")
    c.df <- do.call("rbind", ss.classification)
}


## data format for the classifier:
# rownames should be gene-centric data
# data does not need to be normalized, MS1 quantification is required.
classification <-  ktsp.13(data = data)
