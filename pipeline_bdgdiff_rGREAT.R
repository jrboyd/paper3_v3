source('setup.R')
source("module_bdgdiff_STEIN_ChIPseq.R")
# source("module_ESR1_and_enhancers.R")
source("functions_rGREAT.R")

job1 = wrap_great(subject38 = all_bdgdiff_tables$MCF7_H3K4AC$from_ctrl_to_e2$up_bg,
           query38 = all_bdgdiff_tables$MCF7_H3K4AC$from_ctrl_to_e2$up,
           subject_name = "MCF7_H3K4AC_ctrl_vs_e2",
           query_name = "up_e2",with = T)

job2 = wrap_great(subject38 = all_bdgdiff_tables$MCF7_H3K4AC$from_ctrl_to_e2$down_bg,
                 query38 = all_bdgdiff_tables$MCF7_H3K4AC$from_ctrl_to_e2$down,
                 subject_name = "MCF7_H3K4AC_ctrl_vs_e2",
                 query_name = "down_e2",with = T)

job1 = wrap_great(subject38 = all_bdgdiff_tables$MCF7_H3K4ME3$from_ctrl_to_e2$up_bg,
                  query38 = all_bdgdiff_tables$MCF7_H3K4ME3$from_ctrl_to_e2$up,
                  subject_name = "MCF7_H3K4ME3_ctrl_vs_e2",
                  query_name = "up_e2",with = T)

job2 = wrap_great(subject38 = all_bdgdiff_tables$MCF7_H3K4ME3$from_ctrl_to_e2$down_bg,
                  query38 = all_bdgdiff_tables$MCF7_H3K4ME3$from_ctrl_to_e2$down,
                  subject_name = "MCF7_H3K4ME3_ctrl_vs_e2",
                  query_name = "down_e2",with = T)

job1 = wrap_great(subject38 = all_bdgdiff_tables$MCF7_H3K27ME3$from_ctrl_to_e2$up_bg,
                  query38 = all_bdgdiff_tables$MCF7_H3K27ME3$from_ctrl_to_e2$up,
                  subject_name = "MCF7_H3K27ME3_ctrl_vs_e2",
                  query_name = "up_e2",with = T)

job2 = wrap_great(subject38 = all_bdgdiff_tables$MCF7_H3K27ME3$from_ctrl_to_e2$down_bg,
                  query38 = all_bdgdiff_tables$MCF7_H3K27ME3$from_ctrl_to_e2$down,
                  subject_name = "MCF7_H3K27ME3_ctrl_vs_e2",
                  query_name = "down_e2",with = T)

job1 = wrap_great(subject38 = all_bdgdiff_tables$MCF7_H3K27AC$from_ctrl_to_e2$up_bg,
                  query38 = all_bdgdiff_tables$MCF7_H3K27AC$from_ctrl_to_e2$up,
                  subject_name = "MCF7_H3K27AC_ctrl_vs_e2",
                  query_name = "up_e2",with = T)

job2 = wrap_great(subject38 = all_bdgdiff_tables$MCF7_H3K27AC$from_ctrl_to_e2$down_bg,
                  query38 = all_bdgdiff_tables$MCF7_H3K27AC$from_ctrl_to_e2$down,
                  subject_name = "MCF7_H3K27AC_ctrl_vs_e2",
                  query_name = "down_e2",with = T)