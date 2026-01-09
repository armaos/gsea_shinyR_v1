ora_tab_list %>%
  walk(~removeTab("tabs_plots_ora_dynamic", .x))
ora_tab_list <<- NULL

