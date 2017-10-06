# Title     : comparison between motifs using the 0'th-order models
# Objective : reduce PWM redandency in PEnG!motif output due to motif similarity scores
# Created by: Anja and modified by wanwan
# Created on: 29.09.17

##### functions #####

### read in pwm/bamm file
# todo: here, why do not we use Bamm code for parsing the PWM and bgmodel?
read_pwm <- function(pwm_file, pwm_order, read_order){
    values = scan(file=pwm_file, what = "numeric",strip.white = TRUE,quiet = TRUE, blank.lines.skip = TRUE)
    # Q: isnt pwm_order always 0?
    range = c(1:(pwm_order+1))
    values_per_position = cumsum(4^range)[pwm_order+1]
    pwm_length = length(values)/values_per_position
    offset = 1
    if(read_order > 0){
        offset = cumsum(4^range)[read_order]+1
    }
    # todo: offset is never used in any case
    # which means, only the 0th-order models are applied.
    pwm = c()
    for( i in 0:(pwm_length-1)){
        position = rbind(as.numeric(values[(i*values_per_position+1):(i*values_per_position+4^(read_order+1))]))
        pwm = rbind(pwm,position)
    }
    return(pwm)
}

### read in background model
read_bg <- function(bg_file, read_order, len){
    values = scan(file = bg_file, skip = (2+read_order), nlines=1, quiet=TRUE)
    bg = matrix(rep(values,len),ncol = 4^(read_order+1), byrow=TRUE)
    # Q: here the bg is extended?
    return(bg)
}

### calculate distance between pwm to background
get_dist_bg<- function(pwm,bg,pwm_log_pwm,bg_log_bg){
    avg = (pwm + bg ) / 2
    ### distance due to eq.188, replace p_prime with p_bg
    dist = 0.5 * rowSums(pwm_log_pwm + bg_log_bg - 2* avg * log2(avg))
    return(dist)
}

### calculate minimum distance between two pwms
get_min_dist_p_q <- function(p, plogp, q, qlogq, min_overlap, bg, bglogbg, pad_underhangs){
    min_d = 10e6
    min_W = 0
    min_off_p = 0
    min_off_q = 0

    len_p = dim(p)[1]
    len_q = dim(q)[1]
    print(len_p)
    print(len_q)
    max_ol = min(len_p,len_q)
    # todo: here, should it not be always the case that min_ol==min_overlap??
    min_ol = min(len_p,len_q,min_overlap)
    off_q = 0
    for( off_p in 0:(len_p-min_ol)){
        len_ol = min(max_ol, len_p-off_p)
        # pad "underhangs"
        if(pad_underhangs){
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
        if(pad_underhangs ){
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

### calcaulate distance between 2 pwms
calc_pwm_dist <- function(p,q, off_p,off_q,len_ol, plogp, qlogq){
    ol = c(1:len_ol)
    avg = (p[off_p + ol,] + q[off_q + ol,])/2
    distance = sum(plogp[off_p + ol,] + qlogq[off_q+ol,] - 2 * avg * log2(avg))
    # todo: is it not supposed to be:
    # due to e.q. 188
    # distance = sum(plogp[off_p + ol,] + qlogq[off_q+ol,] - 2 * avg * log2(avg)) * 0.5
    return(distance)
}

### get score for motif - db comparison
get_Scores <- function(pwm,pwm_log_pwm, bg, bg_log_bg, db_folder, db_order, read_order, alpha, min_overlap, p_bg_dist, p_rev_bg_dist, pwm_rev, pwm_rev_log_pwm_rev, pad_underhangs){
    scores = c()
    names = c()
    db_files = list.files( path = db_folder, full.names=TRUE, pattern = db_file_pattern_ending, recursive = TRUE)

    # go through database and compare each to query motif
    for(db_file in db_files){
        pwm_db = read_pwm(db_file, db_order, read_order)
        pwm_db_log_pwm_db = pwm_db * log2(pwm_db)
        if(dim(pwm_db)[1] != dim(bg)[1]){
            bg = matrix(rep(bg[1,],dim(pwm_db)[1]),ncol = 4^(read_order+1), byrow=TRUE)
            bg_log_bg = bg * log2(bg)
        }
        q_bg_dist    = get_dist_bg( pwm_db, bg, pwm_db_log_pwm_db, bg_log_bg)
        p_q_info     = get_min_dist_p_q( pwm, pwm_log_pwm, pwm_db, pwm_db_log_pwm_db, min_overlap, bg, bg_log_bg, pad_underhangs)
        p_rev_q_info = get_min_dist_p_q( pwm_rev, pwm_rev_log_pwm_rev, pwm_db, pwm_db_log_pwm_db, min_overlap, bg, bg_log_bg, pad_underhangs)

        # calc matching score either on original db motif or on reverse complement, whatever has the smaller distance to the query motif
        if(p_q_info["dist"]<p_rev_q_info["dist"]){
            if(pad_underhangs){
                q_bg = sum(q_bg_dist)
                p_bg = sum(p_bg_dist)
            }else{
                q_bg = sum(q_bg_dist[p_q_info["off_q"]:(p_q_info["off_q"]+p_q_info["W"]-1)])
                p_bg = sum(p_bg_dist[p_q_info["off_p"]:(p_q_info["off_p"]+p_q_info["W"]-1)])
            }
            # due to eq. 187
            s_p_q = alpha * (p_bg +q_bg)  - p_q_info["dist"]
            line_in = c(s_p_q, p_q_info["dist"], p_q_info["OL"], p_q_info["off_p"],p_q_info["off_q"], q_bg, p_bg,p_q_info["W"], 0, paste0(unlist(strsplit(basename(db_file),"_"))[1], collapse="_") )
        }else{
            if(pad_underhangs){
                q_bg = sum(q_bg_dist)
                p_bg = sum(p_bg_dist)
            }else{
                q_bg = sum(q_bg_dist[p_rev_q_info["off_q"]:(p_rev_q_info["off_q"]+p_rev_q_info["W"]-1)])
                p_bg = sum(p_bg_dist[p_rev_q_info["off_p"]:(p_rev_q_info["off_p"]+p_rev_q_info["W"]-1)])
            }
            # due to eq. 187
            s_p_q = alpha * (p_bg + q_bg) - p_rev_q_info["dist"]
            line_in = c(s_p_q, p_rev_q_info["dist"],p_rev_q_info["OL"], p_rev_q_info["off_p"], p_rev_q_info["off_q"],q_bg, p_bg, p_rev_q_info["W"], 1, paste0(unlist(strsplit(basename(db_file),"_"))[1], collapse="_"))
        }
        scores = rbind(scores, line_in)
    }
    colnames(scores)<- c("score", "q_p", "OL","off_M1","off_M2", "q_bg", "p_bg", "W","reverseComplement", "Target.ID")
    return(scores)
}

## calculate p- and e-value based on shuffled score distribution
calculate_p_and_e_value <- function(info_real, info_fake, e){

    if(dim(info_real)[1]<2){
        info_real = cbind(info_real,
        "FP"=rep(0,max(1,dim(info_real)[1])),
        "Sl_higher"= sapply(as.numeric(info_real["score"]), FUN = function(x){min(as.numeric(info_fake[as.numeric(info_fake[,"score"])>=x,"score"]))}),
        "Sl_lower" = sapply(as.numeric(info_real["score"]), FUN = function(x){max(as.numeric(info_fake[as.numeric(info_fake[,"score"])<=x,"score"]))}) )

    }else{
        info_real = cbind(info_real,
        "FP"=rep(0,max(1,dim(info_real)[1])),
        "Sl_higher"= sapply(as.numeric(info_real[,"score"]), FUN = function(x){min(as.numeric(info_fake[as.numeric(info_fake[,"score"])>=x,"score"]))}),
        "Sl_lower" = sapply(as.numeric(info_real[,"score"]), FUN = function(x){max(as.numeric(info_fake[as.numeric(info_fake[,"score"])<=x,"score"]))}) )
    }
    info_fake = cbind(info_fake, "FP"=rep(1,max(1,dim(info_fake)[1])),
    "Sl_higher" = info_fake[,"score"],
    "Sl_lower"  = info_fake[,"score"])
    info_fake = info_fake[order(as.numeric(info_fake[,"score"]), decreasing=TRUE),]

    info_all = rbind(info_real,info_fake)
    info_all = info_all[order(as.numeric(info_all[,"score"]), decreasing = TRUE),]
    info_all = cbind(info_all, "FPl"=cumsum(as.numeric(info_all[,"FP"])))

    Nminus = dim(info_fake)[1]
    info_all = cbind(info_all, "p-value"=(as.numeric(info_all[,"FPl"])/Nminus + (1/Nminus) * (as.numeric(info_all[,"Sl_higher"]) - as.numeric(info_all[,"score"] )) / ( as.numeric(info_all[, "Sl_higher"]) - as.numeric(info_all[,"Sl_lower"]) + e)))

    ntop = min(50, (0.1*Nminus))
    Sntop = as.numeric(info_fake[ntop,"score"])
    lambda = (1/ntop) * sum(as.numeric(info_fake[1:ntop,"score"]) - Sntop)
    Sntop_i = max(which(as.numeric(info_all[,"score"])==Sntop))
    for(i in 1:Sntop_i){
        # interpolate for highest p-values in order to avoid false calculations due to less datapoints
        info_all[i,"p-value"] = (ntop/Nminus)*exp(-(as.numeric(info_all[i,"score"]) - Sntop)/lambda)
    }
    for(i in c(which(is.na(info_all[,"Sl_lower"])))){
        # Sl_lower is not defined
        info_all[i,"p-value"] = 1
    }

    info_all = cbind(info_all, "e-value"=as.numeric(info_all[,"p-value"])*dim(info_real)[1])

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

#######################################################
#### BaMM-compare for finding motif-motif matches
#######################################################

# default parameter settings
alpha                = 1
min_overlap          = 4
pwm_order            = 4
db_order             = 4
read_order           = 0
shuffle_times        = 50
p_val_limit          = 0.01
e                    = 1e-5
pad_underhangs       = TRUE
db_folder = '/home/wanwan/benchmark/foo/filterPWM'
query_folder = db_folder
outDir = paste0("/home/wanwan/benchmark/foo/filterPWM")


db_file_pattern_ending         = ".ihbcp"
query_file_pattern_ending      = ".ihbcp"
background_file_pattern_ending = ".hbcp"

args <- commandArgs(trailingOnly = TRUE)
splits <- strsplit(args,split='=')
for(split in splits) {
    if(split[1] == '--outDir'        ){outDir        <- split[2]} # where the output shall be stored
    else if(split[1] == '--dbDir'    ){db_folder     <- split[2]} # folder of db to compare against        ( expected .ihbcp ending "db_file_pattern_ending")
    else if(split[1] == '--queryDir' ){query_folder  <- split[2]} # folder of query motifs to compare with ( expected .ihbcp ending "query_file_pattern_ending")
    else if(split[1] == '--dbOrder'  ){db_order      <- split[2]} # order of the motifs in the db
    else if(split[1] == '--qOrder'   ){pwm_order     <- split[2]} # order of the query motifs
    else if(split[1] == '--readOrder'){read_order    <- split[2]} # order on which to do the comparison ( currently only implemented for 0-th order comparison )
    else if(split[1] == '--p_value'  ){p_val_limit   <- split[2]} # cutoff to report significant motif-motif comparisons
    else if(split[1] == '--sampling' ){shuffle_times <- split[2]} # amount of shuffled negatives per query motif to compute p-values on (20 should be fine, 50 is better but takes long)
}

#outDir          = as.numeric(outDir)
p_val_limit     = as.numeric(p_val_limit)
shuffle_times   = as.numeric(shuffle_times)
read_order      = as.numeric(read_order)
db_order        = as.numeric(db_order)
pwm_order       = as.numeric(pwm_order)

pwm_files = list.files(path = query_folder, full.names = TRUE, pattern = query_file_pattern_ending, recursive = TRUE)
name_list = list.files(path = query_folder, full.names = FALSE, pattern = query_file_pattern_ending, recursive = TRUE)

counter = 0

for (pwm_filepath in pwm_files){
    counter = counter + 1
    message(" PWM #", counter, '\t', name_list[counter])
    # build PWM matrix
    pwm     = read_pwm(pwm_filepath, pwm_order, read_order)
    pwm_log_pwm = pwm *log2(pwm)
    # build PWM matrix for reverse complement
    pwm_rev =  pwm[dim(pwm)[1]:1,dim(pwm)[2]:1]
    pwm_rev_log_pwm_rev = pwm_rev * log2(pwm_rev)

    # build bg matrix
    prefix = paste(unlist(strsplit(pwm_files, "_"))[1], collapse="_")
    bg_filepath = paste0(prefix, background_file_pattern_ending )
    bg      = read_bg(bg_filepath, read_order, dim(pwm)[1])
    bg_log_bg  = bg *log2(bg)
    p_bg_dist = get_dist_bg( pwm, bg, pwm_log_pwm, bg_log_bg)
    # build bg matrix for reverse complement
    p_rev_bg_dist = get_dist_bg( pwm_rev, bg, pwm_rev_log_pwm_rev, bg_log_bg)

    print( pad_underhangs)

    info_real = get_Scores(pwm,pwm_log_pwm, bg, bg_log_bg, db_folder, db_order, read_order, alpha, min_overlap, p_bg_dist, p_rev_bg_dist, pwm_rev, pwm_rev_log_pwm_rev, pad_underhangs)

    # get background scores on shuffled pwms (inclusive reverseComplement)
    Worst_Dist = get_min_dist_p_q(pwm, pwm_log_pwm, bg, bg_log_bg, min_overlap, bg, bg_log_bg, pad_underhangs)

    info_fake = c()
    for( x in 1:shuffle_times ){
        message("Shuffle #", x)

        # shuffle until randomized motif is only to 50% matchable to original motif
        not_shuffled_enough = TRUE
        shuff_counter = 0
        while(not_shuffled_enough){
            shuff_counter    = shuff_counter + 1
            shuff            = shuffle_pwm(pwm)
            shuff_log_shuff  = shuff * log2( shuff )
            shuff_dist = get_min_dist_p_q(pwm, pwm_log_pwm, shuff, shuff_log_shuff, min_overlap, bg, bg_log_bg,pad_underhangs)
            shuff_rev               = shuff[dim(shuff)[1]:1,dim(shuff)[2]:1]
            shuff_rev_log_shuff_rev = shuff_rev * log2(shuff_rev)
            shuff_dist_rev = get_min_dist_p_q(pwm, pwm_log_pwm, shuff_rev, shuff_rev_log_shuff_rev, min_overlap, bg, bg_log_bg,pad_underhangs)

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

        info_shuffle = get_Scores(shuff,shuff_log_shuff,bg, bg_log_bg, db_folder, db_order, read_order, alpha, min_overlap, p_bg_dist, p_rev_bg_dist, shuff_rev, shuff_rev_log_shuff_rev, pad_underhangs)
        info_fake = rbind(info_fake, info_shuffle)
    }

    # calculate p_value and e_value
    info_all = calculate_p_and_e_value(info_real, info_fake , e)
    info_real <- info_all[info_all[,"FP"]== 0,]
    info_real <- cbind(info_real, "Query.ID"=rep(x = basename(name), dim(info_real)[1]))

    best_matches = info_real[which(info_real[,"p-value"]<= p_val_limit),]
    # output results
    if(dim(best_matches)[1] == 0){
        message('no matches!')
    }else{
        for(i in c(1:dim(best_matches)[1])){
            message(rownames(best_matches)[i] , ' ' , best_matches[i,"p-value"] , ' ' , best_matches[i,"e-value"] , ' ' , best_matches[i,"score"] , ' ' , best_matches[i,"W"])
        }
    }

}

