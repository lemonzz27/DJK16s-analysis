Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")

.libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(scales)
})

OUT <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "outputs/20261002_manuscript_revision/RF")

importance_all <- read.csv(file.path(OUT, "RF_importance_All_ConF_USF.csv"), check.names = FALSE)

importance_all$Sig[is.na(importance_all$Sig)] <- ""

importance_all$Group <- factor(importance_all$Group, levels = c("All", "ConF", "USF"))

importance_all$Label <- sprintf("%.2f %s", importance_all$`%IncMSE`, importance_all$Sig)

make_plot <- function(groups, base_size = 9) {
  pd <- importance_all %>%
    filter(Group %in% groups) %>%
    group_by(Group) %>%
    arrange(`%IncMSE`, .by_group = TRUE) %>%
    mutate(Variable_facet = paste(Variable, Group, sep = "___")) %>%
    ungroup()
  pd$Variable_facet <- factor(pd$Variable_facet, levels = pd$Variable_facet)
  max_abs <- max(abs(pd$`%IncMSE`))
  pd$Label_y <- pd$`%IncMSE` + ifelse(pd$`%IncMSE` >= 0, 1, -1) * 0.025 * max_abs
  pd$Label_hjust <- ifelse(pd$`%IncMSE` >= 0, 0, 1)
  ggplot(pd, aes(x = Variable_facet, y = `%IncMSE`, fill = `%IncMSE`)) + geom_hline(yintercept = 0,
    linewidth = 0.35, color = "black") + geom_col(width = 0.78, color = "black", linewidth = 0.25) +
    geom_text(aes(y = Label_y, label = Label, hjust = Label_hjust), size = base_size/ggplot2::.pt,
      color = "black", family = "Arial") + scale_y_continuous(expand = expansion(mult = c(0.15,
    0.22))) + coord_flip(clip = "off") + facet_wrap(~Group, ncol = 1, scales = "free_y") + scale_x_discrete(labels = function(x) sub("___.*$",
    "", x)) + scale_fill_gradientn(colors = c(brewer.pal(5, "Blues")[2], brewer.pal(5, "YlGnBu")[1],
    brewer.pal(6, "OrRd")[4]), values = rescale(c(min(pd$`%IncMSE`), median(pd$`%IncMSE`), max(pd$`%IncMSE`))),
    guide = "none") + labs(x = NULL, y = "Permutation importance (scaled %IncMSE)") + theme_test(base_size = base_size,
    base_family = "Arial") + theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white",
    color = "black", linewidth = 0.4), strip.text = element_text(face = "bold", size = base_size +
    1), axis.text = element_text(color = "black", size = base_size), axis.title.x = element_text(size = base_size +
    1), plot.margin = ggplot2::margin(5, 8, 5, 5))
}

export_plot <- function(p, name, width_mm, height_mm) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = width_mm, height = height_mm, units = "mm",
    dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".svg")), p, width = width_mm, height = height_mm, units = "mm",
    device = svglite::svglite, bg = "white")
}

export_plot(make_plot(c("ConF", "USF")), "Figure_S10_RF_ConF_USF_without_Q", 170, 204)

export_plot(make_plot(c("All", "ConF", "USF")), "RF_importance_All_ConF_USF_without_Q", 183, 275)

for (g in c("All", "ConF", "USF")) export_plot(make_plot(g), paste0("RF_importance_", g, "_without_Q"),
  170, 106)
