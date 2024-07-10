#Commented lines are for figures not being generated currently

# ggsave(paste0(fig_folder, "/le_both.pdf"), le_both, width = width, height = 6.475/5.55*width)
# ggsave(paste0(fig_folder, "/le_optimal.pdf"), le2.plot, width = width, height=3.4/5.55*width)
ggsave(paste0(fig_folder, "/le_baselines.png"),le_baselines, width = 1.3*width, height=1.3*width)
# ggsave(paste0(fig_folder, "/le_continent_mob.png"),le_continent_mob, width = 1.3*width, height=1.3*width)
ggsave(paste0(fig_folder, "/bsa_sens_takeup.png"),df_sens_takeup.plot, width = 1.3*width, height=1.3*width)
# ggsave(paste0(fig_folder, "/burden_transmission_eff.pdf"), burden_vetrans+
#          theme(legend.spacing.x = unit(0.1, 'in'), 
#                text = element_text(size=14),
#                legend.text = element_text(size=10),
#                legend.key.size = unit(0.4, "cm")), width = width, height=0.75*width)
# ggsave(paste0(fig_folder, "/delay_both.pdf"), delay_both, width = width, height = 6.475/5.55*width)
# ggsave(paste0(fig_folder, "/delay_optimal.pdf"), delay_optimal.fig, width = width, height = 3.4/5.55*width)

# ggsave(paste0(fig_folder, "/delay_switch.pdf"), delay_switch + theme(text = element_text(size=9)), width = width, height = 3.7/5.55*width)
# ggsave(paste0(fig_folder, "/g2.pdf"), fig_g2, width = 0.85*width, height=0.85*width)
# ggsave(paste0(fig_folder, "/g2_reductions_only.pdf"), g2b, width = width, height=0.6*width)
# ggsave(paste0(fig_folder, "/sgg_age.pdf"), sgg_age + theme(text = element_text(size=9)), width = width, height=1.85/5.55*width)
ggsave(paste0(fig_folder, "/g1_joint.png"),g1_joint, width = width, height=1.2*width)
ggsave(paste0(fig_folder, "/g2_comp.png"),gg2_comp, width = width, height=1.2*width)
ggsave(paste0(fig_folder, "/g1_comp_joint.png"),g1_comp_joint, width = width, height=1.2*width)
