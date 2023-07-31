pre_proc_data_plot_trim <- function(data,
                               inters_T2_T3) {
  
  data_to_plot <- data %>% 
    mutate(Exposure = ifelse(Exposure == "t2_log2", "Urine pool 1", ifelse(Exposure == "t3_log2", "Urine pool 2", Exposure)),
           sign = case_when(Exposure == "Urine pool 1" & raw_p_value_SEPAGES < 0.05 ~ "<0.05 pool 1", 
                            Exposure == "Urine pool 2" & raw_p_value_SEPAGES < 0.05 ~ "<0.05 pool 2",
                            TRUE ~ ">0.05"),
           sign = ifelse(CpG %in% inters_T2_T3, "<0.05 pool 1 and pool 2", sign),
           sign = factor(sign, levels = c("<0.05 pool 1 and pool 2", "<0.05 pool 1", "<0.05 pool 2", ">0.05")), 
           
           Gene = ifelse(sign != ">0.05", Gene, ""),
           Gene = ifelse(sign == "<0.05 pool 1 and pool 2" & Exposure == "Urine pool 2", "", Gene),
           Gene = stringr::word(Gene, 1, sep = ";"),
           Gene = ifelse(Gene == "Unknown", "Interegenic reg.", Gene))
  
  return(data_to_plot)
}

plot_data <- function(data,
                      title) {
  plot <- data %>% 
    dplyr::arrange(desc(sign)) %>% 
    dplyr::mutate("Exposure assessment time point" = Exposure) %>% 
    ggplot2::ggplot(ggplot2::aes(x = Estimate, y = Estimate_SEPAGES, label = Gene)) +
    ggplot2::geom_point(ggplot2::aes(shape = `Exposure assessment time point`, color = sign)) +
    ggplot2::scale_color_manual(values = c("purple", 
                                           "red",
                                           "blue",
                                           "lightgray"), name = "p-value")+ 
    ggplot2::guides(shape = ggplot2::guide_legend(order = 1),
                    colour = ggplot2::guide_legend(order = 2)) + 
    ggplot2::geom_smooth(ggplot2::aes(group = Exposure), 
                         method = "lm", 
                         se = FALSE) + 
    ggrepel::geom_text_repel(box.padding = 1, max.overlaps = Inf, fontface = "italic") +
    ggplot2::theme(legend.position = "right") + 
    ggplot2::ylab("β for the current study (SEPAGES)") + 
    ggplot2::xlab("β for the previous study (EDEN)") +
    ggplot2::ylim(-0.012, 0.015) +
    ggplot2::xlim(-0.012, 0.015) +
    ggplot2::theme_bw(base_size = 10) + 
    ggplot2::ggtitle(title)
  
  return(plot)
}