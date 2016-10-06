source("setup.R")
source("configure.R")
source("module_DE_res.R")
source("module_ESR1_binding.R")
# source("module_MCF7_enhancers.R")

# GSE40129$MCF7_e2_H3K4ME1
# GSE40129$MCF7_ctrl_H3K4ME1
# MCF7_H3K4ME1_venn = venn_peaks(list(ctrl = GSE40129$MCF7_ctrl_H3K4ME1, 
#                                     e2 = GSE40129$MCF7_e2_H3K4ME1))
# MCF7_H3K4ME1_and_ESR1_venn = venn_peaks(list(ESR1 = GSE40129$MCF7_e2_ESR1, 
#                                              H3K4ME1_e2 = GSE40129$MCF7_e2_H3K4ME1))

e2_olaps = venn_peaks(list(MCF7_e2_H3K4ME1 = GSE40129$MCF7_e2_H3K4ME1, MCF7_e2_H3K27AC = GSE40129$MCF7_e2_H3K27AC))
ctrl_olaps = venn_peaks(list(MCF7_ctrl_H3K4ME1 = GSE40129$MCF7_ctrl_H3K4ME1, MCF7_ctrl_H3K27AC = GSE40129$MCF7_ctrl_H3K27AC))

e2_enh = reduce(c(e2_olaps$`MCF7_e2_H3K4ME1 with MCF7_e2_H3K27AC` , e2_olaps$`MCF7_e2_H3K27AC with MCF7_e2_H3K4ME1`))
ctrl_enh = reduce(c(ctrl_olaps$`MCF7_ctrl_H3K4ME1 with MCF7_ctrl_H3K27AC` , ctrl_olaps$`MCF7_ctrl_H3K27AC with MCF7_ctrl_H3K4ME1`))

enh_change = venn_peaks(list(e2_enh = e2_enh, ctrl_enh = ctrl_enh))
enh_gained = enh_change$`e2_enh no ctrl_enh`
enh_gained_bg = e2_enh
enh_lost = enh_change$`ctrl_enh no e2_enh`
enh_lost_bg = ctrl_enh

enh_gained_job = wrap_great(subject38 = e2_enh, 
                            query38 = ctrl_enh, 
                            subject_name = "e2_enh", 
                            query_name = "ctrl_enh", 
                            with = F)

enh_lost_job = wrap_great(subject38 = ctrl_enh, 
                          query38 = e2_enh, 
                          subject_name = "ctrl_enh", 
                          query_name = "e2_enh", 
                          with = F)

enh_ESR1_bound_job = wrap_great(subject38 = e2_enh, 
                            query38 = GSE40129$MCF7_e2_ESR1, 
                            subject_name = "e2_enh", 
                            query_name = "e2_ESR1", 
                            with = T)

enh_ESR1_absent_job = wrap_great(subject38 = e2_enh, 
                                query38 = GSE40129$MCF7_e2_ESR1, 
                                subject_name = "e2_enh", 
                                query_name = "e2_ESR1", 
                                with = F)

ESR1_at_enh_job = wrap_great(subject38 = GSE40129$MCF7_e2_ESR1, 
                                 query38 = e2_enh, 
                                 subject_name = "ESR1", 
                                 query_name = "e2_enh", 
                                 with = T)

ESR1_not_at_enh_job = wrap_great(subject38 = GSE40129$MCF7_e2_ESR1, 
                             query38 = e2_enh, 
                             subject_name = "ESR1", 
                             query_name = "e2_enh", 
                             with = F)

par(mfrow = c(1, 3))
res_gained = plotRegionGeneAssociationGraphs(enh_gained_job, )
res_lost = plotRegionGeneAssociationGraphs(enh_lost_job)
