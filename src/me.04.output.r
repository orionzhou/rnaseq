#{{{ header
source("me.fun.r")
t_cfg
#}}}

#{{{ write output
diro = file.path(dird, '15_output')
for(i in 1:nrow(t_cfg)) {
    sid = t_cfg$sid[i]; study = t_cfg$study[i]; meta = t_cfg$meta[i]
    #
    diri = file.path(dird, '11_qc', sid)
    fi = file.path(diri, '20.rc.norm.rda')
    if(!file.exists(fi)) next
    cat(sid, study, "\n")
    x = load(fi)
    if(is.na(meta) | meta != T) {
        th = get_read_list(dird, sid)
        tm = tm %>% filter(SampleID %in% th$SampleID)
        res = merge_reps(th, tm, sid)
        th = res$th; tm = res$tm
    }
    #stopifnot(nrow(th) * ngene == nrow(tm))
    #
    fo = sprintf("%s/%s.rda", diro, sid)
    save(th, tm, file = fo)
}
#}}}

