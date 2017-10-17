#!/usr/bin/env Rscript

##### functions #####
options( warn = -1 )
library( rtfbs )        # Transcription Factor Binding Site Identification Tool
library( argparse )     # parse arguments

### read in meme formated files using rtfbs package
read_meme <- function(pwm_file){
    pwm_length = 0
    # read in and log-transform
    pwms = read.pwm(pwm_file)
    # check whether one or more pwms are in the file
    if( class(pwms) == "list" ){
        # exp to remove the log transformation
        for(i in c(1:length(pwms))){
            pwms[[i]] <- exp(pwms[[i]])
            if( dim(pwms[[i]])[1] > pwm_length){
                pwm_length = dim(pwms[[i]])[1]
            }
        }
    } else {
        pwm <- exp(pwms)
        pwm_length = dim(pwms)[1]
        pwms <- list()
        pwms[[1]]<-pwm
    }
    # add motif iupac names and extract background frequencies
    lines <- readLines(pwm_file)
    IUPACS <- lines[grep("^MOTIF*", lines)]
    names(pwms) <- unlist(lapply(IUPACS, function(x){unlist(strsplit(x, "\\ "))[2]}))

    Background <- lines[grep("^Background*", lines)+1]
    bg <- matrix(rep(as.numeric(unlist(strsplit(Background, " "))[c(2,4,6,8)])),pwm_length,ncol=4^(0+1), byrow=TRUE)

    return(list(pwms, bg))
}

### read in bamm file
read_pwm <- function(pwm_file, pwm_order, read_order){
    values = scan(file=pwm_file, what = "numeric",strip.white = TRUE,quiet = TRUE, blank.lines.skip = TRUE)
    range = c(1:(pwm_order+1))
    values_per_position = cumsum(4^range)[pwm_order+1]
    pwm_length = length(values)/values_per_position
    offset = 1
    if(read_order > 0){
        offset = cumsum(4^range)[read_order]+1
    }
    pwm = c()
    for( i in 0:(pwm_length-1)){
        position = rbind(as.numeric(values[(i*values_per_position+1):(i*values_per_position+4^(read_order+1))]))
        pwm = rbind(pwm,position)
    }
    pwms <- list()
    pwms[[1]]<-pwm
    names(pwms)[1] <- paste0(unlist(strsplit(basename(pwm_file),"_"))[1], collapse="_")
    return(pwms)
}

### read in background model
read_bg <- function(bg_file, read_order, len){
    values = scan(file = bg_file, skip = (2+read_order), nlines=1, quiet=TRUE)
    bg = matrix(rep(values,len),ncol = 4^(read_order+1), byrow=TRUE)
    return(bg)
}

### calculate distance between pwm to background (eq.188)
get_dist_bg<- function(pwm,bg,pwm_log_pwm,bg_log_bg){

    avg = (pwm + bg[c(1:dim(pwm)[1]),] ) / 2
    dist = 0.5 * rowSums(pwm_log_pwm + bg_log_bg[c(1:dim(pwm)[1]),] - 2* avg * log2(avg))
    return(dist)
}

### calculate minimum distance between two pwms
get_min_dist_p_q <- function(p, plogp, q, qlogq, min_overlap, bg, bglogbg,pad_underhangs){

    min_d = 10e6
    min_W = 0
    min_off_p = 0
    min_off_q = 0

    len_p = dim(p)[1]
    len_q = dim(q)[1]
    max_ol = min(len_p,len_q)
    min_ol = min(len_p,len_q,min_overlap)
    off_q = 0
    for( off_p in 0:(len_p-min_ol)){
        len_ol = min(max_ol, len_p-off_p)
        # pad "underhangs"
        if( pad_underhangs ){
            max_len = len_q+len_p-len_ol

            # extend p  = off_q_BG + p + (max_len-len_p-off_q)_BG
            p_padded = rbind(matrix(rep(bg[1,],off_q),ncol = 4^(read_order+1), byrow=TRUE), p,
                            matrix(rep(bg[1,],(max_len-len_p-off_q)),ncol = 4^(read_order+1), byrow=TRUE))

            plogp_padded = rbind(matrix(rep(bglogbg[1,],off_q),ncol = 4^(read_order+1), byrow=TRUE), plogp,
                                matrix(rep(bglogbg[1,],(max_len-len_p-off_q)),ncol = 4^(read_order+1), byrow=TRUE))

            # extend q  = off_p_BG + q + (max_len-len_q-off_p)_BG
            q_padded = rbind(matrix(rep(bg[1,],off_p),ncol = 4^(read_order+1), byrow=TRUE), q,
                            matrix(rep(bg[1,],(max_len-len_q-off_p)),ncol = 4^(read_order+1), byrow=TRUE))

            qlogq_padded = rbind(matrix(rep(bglogbg[1,],off_p),ncol = 4^(read_order+1), byrow=TRUE), qlogq,
                            matrix(rep(bglogbg[1,],(max_len-len_q-off_p)),ncol = 4^(read_order+1), byrow=TRUE))

            distance = calc_pwm_dist(p_padded,q_padded,0,0,max_len,plogp_padded,qlogq_padded)
        } else {
            max_len = len_ol
            distance = calc_pwm_dist(p,q,off_p,0,len_ol, plogp,qlogq)
        }

        if( distance < min_d ){
            min_d = distance
            min_W = len_ol
            min_off_p = off_p
            maxLen = max_len
        }
    }
    off_p = 0
    for( off_q in 0:(len_q-min_ol)){
        len_ol = min(max_ol, len_q-off_q)
        # pad "underhangs"
        if( pad_underhangs ){
            max_len = len_q+len_p-len_ol

            # extend p  = off_q_BG + p + (max_len-len_p-off_q)_BG
            p_padded = rbind(matrix(rep(bg[1,],off_q),ncol = 4^(read_order+1), byrow=TRUE), p,
                            matrix(rep(bg[1,],(max_len-len_p-off_q)),ncol = 4^(read_order+1), byrow=TRUE))

            plogp_padded = rbind(matrix(rep(bglogbg[1,],off_q),ncol = 4^(read_order+1), byrow=TRUE),plogp,
                                matrix(rep(bglogbg[1,],(max_len-len_p-off_q)),ncol = 4^(read_order+1), byrow=TRUE))

            # extend q  = off_p_BG + q + (max_len-len_q-off_p)_BG
            q_padded = rbind(matrix(rep(bg[1,],off_p),ncol = 4^(read_order+1), byrow=TRUE), q,
                            matrix(rep(bg[1,],(max_len-len_q-off_p)),ncol = 4^(read_order+1), byrow=TRUE))

            qlogq_padded = rbind(matrix(rep(bglogbg[1,],off_p),ncol = 4^(read_order+1), byrow=TRUE), qlogq,
                                matrix(rep(bglogbg[1,],(max_len-len_q-off_p)),ncol = 4^(read_order+1), byrow=TRUE))

            distance = calc_pwm_dist(p_padded,q_padded,0,0,max_len,plogp_padded,qlogq_padded)
        }else{
            max_len = len_ol
            distance = calc_pwm_dist(p,q,0,off_q,len_ol,plogp, qlogq)
        }

        if( distance < min_d ){
            min_d = distance
            min_W = len_ol
            min_off_p = 0
            min_off_q = off_q
            maxLen = max_len
        }
    }

    info = c(min_d, min_W, min_off_p, min_off_q, maxLen)
    names(info)<- c("dist", "OL", "off_p","off_q", "W")
    return(info)
}

### calcaulate distance between 2 pwms with the same length
calc_pwm_dist <- function(p,q, off_p,off_q,len_ol, plogp, qlogq){
    ol = c(1:len_ol)
    avg = (p[off_p + ol,] + q[off_q + ol,])/2
    distance = sum(plogp[off_p + ol,] + qlogq[off_q+ol,] - 2 * avg * log2(avg)) * 0.5
    return(distance)
}

### read pwms from db folder
read_db_folder <- function(db_folder, db_order, read_order){
    db_files = list.files( path = db_folder, full.names=TRUE, pattern = db_file_pattern_ending, recursive = TRUE)
    db_pwms = list()
    for(i in c(1:length(db_files))){
        db_pwms[i] = read_pwm(db_files[i], db_order, read_order)
        names(db_pwms)[i] <- paste0(unlist(strsplit(basename(db_files[i]),"_"))[1], collapse="_")
    }
    return(db_pwms)
}

### get score for motif - db comparison
get_Scores <- function(pwm,pwm_log_pwm, bg, bg_log_bg,
                        db_motifs, read_order, alpha, min_overlap,
                        p_bg_dist, p_rev_bg_dist, pwm_rev, pwm_rev_log_pwm_rev, pad_underhangs){
    scores = c()

    for( db_id in names(db_motifs)){
        pwm_db = db_motifs[[db_id]]
        # add little pseudos to pwm_db in case the counts are zero to avoid log problems
        pwm_db[pwm_db == 0]<-1e-5

        pwm_db_log_pwm_db = pwm_db * log2(pwm_db)
        if(dim(pwm_db)[1] != dim(bg)[1]){
            bg = matrix(rep(bg[1,],dim(pwm_db)[1]),ncol = 4^(read_order+1), byrow=TRUE)
            bg_log_bg = bg * log2(bg)
        }
        q_bg_dist     = get_dist_bg( pwm_db, bg, pwm_db_log_pwm_db, bg_log_bg)
        p_q_info      = get_min_dist_p_q( pwm, pwm_log_pwm, pwm_db, pwm_db_log_pwm_db,
                                            min_overlap, bg, bg_log_bg, pad_underhangs)
        p_rev_q_info  = get_min_dist_p_q( pwm_rev, pwm_rev_log_pwm_rev, pwm_db, pwm_db_log_pwm_db,
                                            min_overlap, bg, bg_log_bg, pad_underhangs)

        # calc matching score eitehr on original db motif or on reverse complement, whatever has the smaler distance to the query motif
        if(p_q_info["dist"]<p_rev_q_info["dist"]){
            if(pad_underhangs){
                q_bg = sum(q_bg_dist)
                p_bg = sum(p_bg_dist)
            }else{
                q_bg = sum(q_bg_dist[p_q_info["off_q"]:(p_q_info["off_q"]+p_q_info["W"]-1)])
                p_bg = sum(p_bg_dist[p_q_info["off_p"]:(p_q_info["off_p"]+p_q_info["W"]-1)])
            }
            s_p_q = alpha * (p_bg +q_bg)  - p_q_info["dist"]
            line_in = c(s_p_q, p_q_info["dist"], p_q_info["OL"], p_q_info["off_p"],p_q_info["off_q"],
                        q_bg, p_bg,p_q_info["W"], 0, db_id )
        }else{
            if(pad_underhangs){
                q_bg = sum(q_bg_dist)
                p_bg = sum(p_rev_bg_dist)
            }else{
                q_bg = sum(q_bg_dist[p_rev_q_info["off_q"]:(p_rev_q_info["off_q"]+p_rev_q_info["W"]-1)])
                p_bg = sum(p_bg_dist[p_rev_q_info["off_p"]:(p_rev_q_info["off_p"]+p_rev_q_info["W"]-1)])
            }
            s_p_q = alpha * (p_bg + q_bg) - p_rev_q_info["dist"]
            line_in = c(s_p_q, p_rev_q_info["dist"],p_rev_q_info["OL"], p_rev_q_info["off_p"], p_rev_q_info["off_q"],
                        q_bg, p_bg, p_rev_q_info["W"], 1, db_id)
        }
        scores = rbind(scores, line_in)
    }
    if(class(scores) == "matrix"){
        colnames(scores)<- c("score", "q_p", "OL","off_M1","off_M2", "q_bg", "p_bg", "W","reverseComplement", "Target.ID")
    }else{
        names(scores)<- c("score", "q_p", "OL","off_M1","off_M2", "q_bg", "p_bg", "W","reverseComplement", "Target.ID")
    }

    return(scores)
}

## calculate p- and e-value based on shuffled score distribution
calculate_p_and_e_value <- function(info_real, info_fake, e){

    if(dim(info_real)[1]<2){
        info_real = cbind(
            info_real,
            "FP"=rep(0,max(1,dim(info_real)[1])),
            "Sl_higher"= sapply(as.numeric(info_real["score"]),
                            FUN = function(x){min(as.numeric(info_fake[as.numeric(info_fake[,"score"])>=x,"score"]))}),
            "Sl_lower" = sapply(as.numeric(info_real["score"]),
                            FUN = function(x){max(as.numeric(info_fake[as.numeric(info_fake[,"score"])<=x,"score"]))})
        )

    }else{
        info_real = cbind(
            info_real,
            "FP"=rep(0,max(1,dim(info_real)[1])),
            "Sl_higher"= sapply(as.numeric(info_real[,"score"]),
                            FUN = function(x){min(as.numeric(info_fake[as.numeric(info_fake[,"score"])>=x,"score"]))}),
            "Sl_lower" = sapply(as.numeric(info_real[,"score"]),
                            FUN = function(x){max(as.numeric(info_fake[as.numeric(info_fake[,"score"])<=x,"score"]))})
        )
    }
    info_fake = cbind(
        info_fake,
        "FP"=rep(1,max(1,dim(info_fake)[1])),
        "Sl_higher" = info_fake[,"score"],
        "Sl_lower"  = info_fake[,"score"]
    )
    info_fake = info_fake[order(as.numeric(info_fake[,"score"]), decreasing=TRUE),]

    info_all = rbind(info_real,info_fake)
    info_all = info_all[order(as.numeric(info_all[,"score"]), decreasing = TRUE),]
    info_all = cbind(info_all, "FPl"=cumsum(as.numeric(info_all[,"FP"])))

    Nminus = dim(info_fake)[1]
    info_all = cbind(
                info_all,
                # calculate p-value according to eq. 182
                "p-value"= ( as.numeric(info_all[,"FPl"])/Nminus +
                            (1/Nminus) * (as.numeric(info_all[,"Sl_higher"]) - as.numeric(info_all[,"score"] )) /
                            ( as.numeric(info_all[, "Sl_higher"]) - as.numeric(info_all[,"Sl_lower"]) + e) )
    )

    ntop = min(100, (0.1*Nminus))
    Sntop = as.numeric(info_fake[ntop,"score"])
    lambda = (1/ntop) * sum(as.numeric(info_fake[1:ntop,"score"]) - Sntop)
    Sntop_i = max(which(as.numeric(info_all[,"score"])==Sntop))
    for(i in 1:Sntop_i){
        # interpolate for highest p-values in order to avoid false calculations due to less datapoints
        info_all[i,"p-value"] = ( ntop/Nminus ) * exp(-(as.numeric(info_all[i,"score"]) - Sntop)/lambda)
    }
    for(i in c(which(is.na(info_all[,"Sl_lower"])))){
        # Sl_lower is not defined
        info_all[i,"p-value"] = 1
    }

    info_all = cbind(info_all, "e-value"=as.numeric(info_all[,"p-value"]) * dim(info_real)[1])

    return(info_all)

}

# shuffle pwm in a localized fashion
shuffle_pwm <- function(pwm){
    for(i in c(1:dim(pwm)[1])){
        # 1.) switch A->T with p(0.5) per column
        if(rnorm(1,0,1)>0){
            A = pwm[i,1]
            pwm[i,1]<-pwm[i,4]
            pwm[i,4]<- A
        }
        # 2.) switch C->G with p(0.5) per column
        if(rnorm(1,0,1)>0){
            C = pwm[i,2]
            pwm[i,2]<-pwm[i,3]
            pwm[i,3]<- C
        }
    }
    # locally shuffle positions
    j = c(1:dim(pwm)[1])
    Rj = -j +2*rnorm(dim(pwm)[1],0,1)
    return(pwm[order(Rj, decreasing = TRUE),])
}

# run MMCompare for one Query Motif
MMcompare <- function(pwm,query_name, bg, db_motifs,
                        read_order, shuffle_times, pad_underhangs,
                        p_val_limit, min_overlap, alpha, e, outFile){
    # add little pseudos to pwm in case the counts are zero to avoid log problems
    pwm[pwm == 0]<-1e-5

    pwm_log_pwm         = pwm * log2(pwm)
    bg_log_bg           = bg * log2(bg)
    p_bg_dist           = get_dist_bg( pwm, bg, pwm_log_pwm, bg_log_bg)
    pwm_rev             = pwm[dim(pwm)[1]:1,dim(pwm)[2]:1]
    pwm_rev_log_pwm_rev = pwm_rev * log2(pwm_rev)
    p_rev_bg_dist       = get_dist_bg( pwm_rev, bg, pwm_rev_log_pwm_rev, bg_log_bg)

    info_real = get_Scores(pwm,pwm_log_pwm, bg, bg_log_bg,
                            db_motifs, read_order, alpha, min_overlap,
                            p_bg_dist, p_rev_bg_dist, pwm_rev, pwm_rev_log_pwm_rev, pad_underhangs)

    # get background scores on shuffled pwms (inclusive reverseComplement)
    Worst_Dist = get_min_dist_p_q(pwm, pwm_log_pwm, bg, bg_log_bg, min_overlap, bg, bg_log_bg, pad_underhangs)

    info_fake = c()
    for( x in 1:shuffle_times ){
        #message("Shuffle #", x)

        # shuffle until randomized motif is only to 50% matchable to original motif
        not_shuffled_enough = TRUE
        shuff_counter = 0
        while(not_shuffled_enough){
            shuff_counter           = shuff_counter + 1
            shuff                   = shuffle_pwm(pwm)
            shuff_log_shuff         = shuff * log2( shuff )
            shuff_dist              = get_min_dist_p_q(pwm, pwm_log_pwm, shuff, shuff_log_shuff,
                                                        min_overlap, bg, bg_log_bg, pad_underhangs)
            shuff_rev               = shuff[dim(shuff)[1]:1,dim(shuff)[2]:1]
            shuff_rev_log_shuff_rev = shuff_rev * log2(shuff_rev)
            shuff_dist_rev          = get_min_dist_p_q(pwm, pwm_log_pwm, shuff_rev, shuff_rev_log_shuff_rev,
                                                        min_overlap, bg, bg_log_bg, pad_underhangs)

            if(min(shuff_dist[1],shuff_dist_rev[1]) > Worst_Dist[1]){
                not_shuffled_enough = FALSE
            }
            if(shuff_counter == 500){
                not_shuffled_enough = FALSE
            }
        }
        shuff_log_shuff         = shuff * log2( shuff )
        shuff_rev               = shuff[dim(shuff)[1]:1,dim(shuff)[2]:1]
        shuff_rev_log_shuff_rev = shuff_rev * log2(shuff_rev)

        info_shuffle = get_Scores(shuff, shuff_log_shuff, bg, bg_log_bg,
                                    db_motifs, read_order, alpha, min_overlap,
                                    p_bg_dist, p_rev_bg_dist, shuff_rev, shuff_rev_log_shuff_rev, pad_underhangs)
        info_fake = rbind(info_fake, info_shuffle)
    }

    # calculate p_value and e_value
    info_all = calculate_p_and_e_value(info_real, info_fake , e)
    info_real <- info_all[info_all[,"FP"]== 0,]

    if(class(info_real) == "matrix"){
        best_matches = info_real[which(as.numeric(info_real[,"p-value"]) <= p_val_limit),]
    }else{
        best_matches = info_real[which(as.numeric(info_real["p-value"]) <= p_val_limit)]
    }


    if(outFile != "NA"){
        # create if file doesn't exist yet
        if(!file.exists(outFile)){
            write(x='', file=outFile)
        }
        # output results
        if(class(best_matches) == "character"){
            if(length(best_matches) == 0){
                write(x='no matches!\n', file=outFile, append=TRUE)
            }else{
                write(x='Query Target p-value e-value score overlap-length',file = outFile, append = TRUE)
                write(x=paste(query_name,'\t',
                                best_matches["Target.ID"], '\t',
                                best_matches["p-value"], '\t',
                                best_matches["e-value"], '\t',
                                best_matches["score"], '\t',
                                best_matches["W"]),
                                file=outFile, append=TRUE)
                write(x="\n", file=outFile, append=TRUE)
            }
        }
        else{
            if(dim(best_matches)[1] == 0){
                write(x='no matches!\n',file=outFile, append=TRUE)
            }else{
                write(x="Query Target p-value e-value score overlap-length", file=outFile, append=TRUE)
                for(i in c(1:dim(best_matches)[1])){
                    write(x=paste(query_name,'\t',
                                    best_matches[i,"Target.ID"], '\t',
                                    best_matches[i,"p-value"], '\t',
                                    best_matches[i,"e-value"], '\t',
                                    best_matches[i,"score"], '\t',
                                    best_matches[i,"W"]),
                                    file=outFile, append = TRUE )
                }
                write(x="\n", file=outFile, append=TRUE)
            }
        }
    }else{
        # output results
        if(class(best_matches) == "character"){
            if(length(best_matches) == 0){
                message('no matches!\n')
            }else{
                message("Query Target p-value e-value score overlap-length")
                message(query_name,' ',
                        best_matches["Target.ID"], ' ',
                        best_matches["p-value"], ' ',
                        best_matches["e-value"], ' ',
                        best_matches["score"], ' ',
                        best_matches["W"])
                message("\n")
            }
        }
        else{
            if(dim(best_matches)[1] == 0){
                message('no matches!\n')
            }else{
                message("Query Target p-value e-value score overlap-length")
                for(i in c(1:dim(best_matches)[1])){
                    message(query_name,' ',
                            best_matches[i,"Target.ID"], ' ',
                            best_matches[i,"p-value"], ' ',
                            best_matches[i,"e-value"], ' ',
                            best_matches[i,"score"], ' ',
                            best_matches[i,"W"])
                }
                message("\n")
            }
        }
    }
}

#######################################################
#### BaMM-compare for finding motif-motif matches
#######################################################
db_file_pattern_ending=".ihbcp"

db_folder      <- "NA"
db_order       <- 4
pwm_order      <- 4
read_order     <- 0
querybg        <- "NA"
shuffle_times  <- 100
pad_underhangs <- TRUE
p_val_limit    <- 0.01
min_overlap    <- 4
alpha          <- 1
e              <- 1e-5

parser <- ArgumentParser(description="plot the motif distribution")
# positional arguments
parser$add_argument('query', help="full path to query file (either meme (.meme) or bamm (.ihbcp) format)")

# optional arguments needed only for bamm file format (all of those can be neglected when comparing meme pwm list with itself)
parser$add_argument("--dbDir",       type="character", default="NA", help="full path to database, if not given, query will be screed against itself" )
parser$add_argument("--dbOrder",     type="integer",   default=0,    help="order of motifs stored in the database" )
parser$add_argument("--qOrder",      type="integer",   default=0,    help="order of query motifs for bamm format" )
parser$add_argument("--readOrder",   type="integer",   default=0,    help="order on which to compare query and DB (currently only 0-th order implemented!)" )
parser$add_argument("--bg",          type="character", default="NA", help="full path to query background file for bamm format" )
parser$add_argument("--outFile",     type="character", default="NA", help="full path to output file. If none given results are written to the terminal.")

# optional arguments for parameter settings
parser$add_argument("--sampling",    type="integer",    default=20,     help="number of sampled query motifs to calulate p-value" )
parser$add_argument("--padding",     type="logical",    default=TRUE,   help="boolean, if underhangs due to shifting should be padded" )
parser$add_argument("--pValue",      type="double",     default=0.01,   help="limit for reporting motif-motif comparisons" )
parser$add_argument("--minOverlap",  type="integer",    default=4,      help="minimum overlap between query and db motif" )
parser$add_argument("--alpha",       type="double",     default=1,      help="weighting of the motif strength in the overall score" )
parser$add_argument("--epsilon",     type="double",     default=1e-5,   help="small factor to avoid division by 0" )

args           <- parser$parse_args()
query          <- args$query
db_folder      <- args$dbDir
db_order       <- args$dbOrder
pwm_order      <- args$qOrder
read_order     <- args$readOrder
querybg        <- args$bg
outFile        <- args$outFile
shuffle_times  <- args$sampling
pad_underhangs <- args$padding
p_val_limit    <- args$pValue
min_overlap    <- args$minOverlap
alpha          <- args$alpha
e              <- args$epsilon

if( unlist(strsplit(query, "\\."))[-1] == "ihbcp" ){
    pwm         = read_pwm(query, pwm_order, read_order)
    bg          = read_bg(querybg, read_order, dim(pwm[[1]])[1])
}

if( unlist(strsplit(query, "\\."))[-1] == "meme" ){
    info        = read_meme(query)
    pwm         = info[[1]]
    bg          = info[[2]]
}

if(db_folder != "NA"){
    # THIS NEEDS EXTENSION FOR WORKING WITH A DB OF PWMS! (currently only working with the .ihbcp db)
    db_motifs = read_db_folder(db_folder, db_order, read_order)
}else{
    db_motifs = pwm
}

for(name in names(pwm)){
    MMcompare(pwm[[name]],name, bg, db_motifs,
                read_order, shuffle_times, pad_underhangs,
                p_val_limit, min_overlap, alpha, e, outFile)
}
