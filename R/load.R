#' Semiparametric Causal Mediation Modeling of Semi-Competing Risks
#'
#' This function analyzes semicompeting risks data and gives the estimators of direct and indirect effects, along with their variances.
#'
#' @param df a data frame with designated column order: \code{T1}, \code{T2}, \code{d1}, \code{d2}, \code{Z}, \code{X} where X maybe empty or a matrix.
#' @param effect the causal effect to be estimated. Choices are \code{"DE"} or \code{"IE"}. Default is \code{c("DE", "IE")}
#' @param intervention numeric vector with two elements. Default is \code{c(1, 0)}.
#' @param cal_level the level \code{X} should be evaluated at. It can be a numeric vector with \code{length(cal_level) = dim(X)[2]} or \code{"median"}, \code{"mean"}, etc. Default is \code{"median"}.
#' @param myunit the length to be considered as one unit. Default is \code{"raw"}.
#' @param downsample indicates how many consecutive event should be considered as one event. Default is \code{1}.
#' @param get_variance the method to compute the variance. Choices are \code{"a"} or \code{"b"}. Default is \code{"a"}.
#' @param boot_times the times of bootstrap. Default is \code{1000}.
#' @param timer will show the progress. Default is \code{TRUE}.
#' @param sen_ana doing sensitivity analysis or not. Default is \code{FALSE}.
#' @param num_of_cores the number of cores assigned. Default is \code{1}.
#' @param plot_result will show some primary result. Default is \code{FALSE}.
#' @param threshold specifies if a logistic regression converge or not. Default is \code{1e-10}.
#' @return \code{CASCR} returns a list with components specified by \code{effect}.
#' @export
CASCR = function(df, effect = c('DE', 'IE'), intervention = c(1, 0), cal_level = 'median', myunit = 'raw', downsample = 1, sen_ana = FALSE, get_variance = c('asymptotic'), boot_times = 1000, faster_bootstrap = 1, timer = TRUE, num_of_cores = 1, plot_result = FALSE, variance_method = 'new', threshold = 1e-10, HO = FALSE){
  # effect = c('DE', 'IE'); intervention = c(1, 0); cal_level = 'median'; myunit = 'raw'; downsample = 1; sen_ana = FALSE; get_variance = c('asymptotic'); boot_times = 1000; faster_bootstrap = 1; timer = TRUE; num_of_cores = 10; plot_result = FALSE; variance_method = 'new'; threshold = 1e-10
  # protect the original data
  df = df[, c(1, 3, 2, 4, 5:dim(df)[2])]
  dff = df

  #
  df = data_preprocess(dff, myunit, downsample)
  unique_T2 = df$unique_T2
  b0_time = df$b0_time
  b1_time = df$b1_time

  df = df_shift_to_cal_level(df$df, cal_level)
  ana_cal_level = df$cal_level
  df = df$df

  result = estimate_effect(df, effect, intervention, cal_level = ana_cal_level, sen_ana, get_variance, boot_times, timer, num_of_cores, unique_T2, b0_time, b1_time, variance_method, threshold, HO = HO)
  result$cal_level = ana_cal_level

  BootVariance = sum(c('b', 'B', 'boot', 'bootstrap', 'bootstrapping', 'Boot', 'Bootstrap', 'Bootstrapping') %in% get_variance) > 0
  if(BootVariance){
    if(faster_bootstrap > 1){
      faster_bootstrap = ceiling(faster_bootstrap)
      if(myunit != "raw"){
        myunit = myunit * faster_bootstrap
      }
      downsample = downsample * faster_bootstrap
    }
    get_DE = sum(c('d', 'D', 'de', 'De', 'DE', 'direct effect', 'Direct effect', 'Direct Effect') %in% effect) > 0
    get_IE = sum(c('i', 'I', 'ie', 'Ie', 'IE', 'indirect effect', 'Indirect effect', 'Indirect Effect') %in% effect) > 0

    my_eva_time = unique_T2
    if(get_DE){
      boot_DE_mat = matrix(0, boot_times, length(my_eva_time))
      # Q_stat_DE = rep(0, boot_times)
    }
    if(get_IE){
      boot_IE_mat = matrix(0, boot_times, length(my_eva_time))
      # Q_stat_IE = rep(0, boot_times)
    }

    if(num_of_cores > 1){
      require(foreach)

      ## fetch basic parameter
      m = dim(df)[1]

      ## num_of_cores set-up
      cores = num_of_cores
      cl = snow::makeCluster(cores[1])
      my_functions = c("downsample_func", "data_preprocess", "df_shift_to_cal_level", "do_sen_ana", "estimate_alpha", "estimate_effect", "form_matrix", "get_alpha_variance", "get_beta_variance", "get_counterfactual_hazard", "get_pd", "get_position", "compute_variance", "inv_coxinformation", "make_small", "my_basehaz", "my_eva_fun", "my_sort_mat", "mycoxph", "my_rep_row")
      snow::clusterExport(cl, my_functions)
      doSNOW::registerDoSNOW(cl)
      pb = txtProgressBar(max = boot_times, style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      opts = list(progress = progress)

      boot_effect = foreach(i = 1:boot_times, .options.snow = opts, .combine = 'c', .export = my_functions) %dopar%{
        set.seed(2020 + i)
        boot_index = sample(1:m, m, replace = TRUE)
        boot_df = dff[boot_index, ]

        boot_df = data_preprocess(boot_df, myunit, downsample)
        unique_T2 = boot_df$unique_T2
        b0_time = boot_df$b0_time
        b1_time = boot_df$b1_time

        boot_df = df_shift_to_cal_level(boot_df$df, cal_level)
        boot_cal_level = boot_df$cal_level
        boot_df = boot_df$df

        boot_effect = list(estimate_effect(boot_df, effect, intervention, cal_level = boot_cal_level, sen_ana = FALSE, get_variance = NULL, boot_times = 0, timer = FALSE, num_of_cores = FALSE, unique_T2, b0_time, b1_time, variance_method, threshold = 1e-15))
        gc()
        return(boot_effect)
      }
      snow::stopCluster(cl)
      pracma::fprintf('\n')
      for(i in 1:boot_times){
        if(get_DE){
          # Q_stat_DE[i] = boot_effect[[i]]$DE$Q_stat
          boot_DE_mat[i, ] = my_eva_fun(list(boot_effect[[i]]$DE$effect, boot_effect[[i]]$DE$time), my_eva_time, method = "linear")
        }
        if(get_IE){
          # Q_stat_IE[i] = boot_effect[[i]]$IE$Q_stat
          boot_IE_mat[i, ] = my_eva_fun(list(boot_effect[[i]]$IE$effect, boot_effect[[i]]$IE$time), my_eva_time, method = "linear")
        }
      }
    }else{
      ## fetch basic parameter
      m = dim(df)[1]

      if(timer){
        space = 100
        pracma::fprintf('| bootstrap        20        30        40        50        60        70        80        90    100 |\n')
        loop_count = 1:boot_times
        counter_total = boot_times
        cum_bar_num = my_eva_fun(list(1:space, 1:space / space * counter_total), loop_count, rule = '0')
        bar_num = diff(c(0, cum_bar_num))
      }
      i = 1
      for(i in 1:boot_times){
        # print(i)
        set.seed(2020 + i)
        boot_index = sample(1:m, m, replace = TRUE)
        boot_df = dff[boot_index, ]

        boot_df = data_preprocess(boot_df, myunit, downsample)
        unique_T2 = boot_df$unique_T2
        b0_time = boot_df$b0_time
        b1_time = boot_df$b1_time

        boot_df = df_shift_to_cal_level(boot_df$df, cal_level)
        boot_cal_level = boot_df$cal_level
        boot_df = boot_df$df

        boot_effect = estimate_effect(boot_df, effect, intervention, boot_cal_level, sen_ana = FALSE, get_variance = NULL, boot_times = 0, timer = FALSE, num_of_cores = FALSE, unique_T2, b0_time, b1_time, variance_method, threshold = 1e-15)
        if(get_DE){
          # Q_stat_DE[i] = boot_effect$DE$Q_stat
          boot_DE_mat[i, ] = my_eva_fun(list(boot_effect$DE$effect, boot_effect$DE$time), my_eva_time, method = "linear")
        }
        if(get_IE){
          # Q_stat_IE[i] = boot_effect$IE$Q_stat
          boot_IE_mat[i, ] = my_eva_fun(list(boot_effect$IE$effect, boot_effect$IE$time), my_eva_time, method = "linear")
        }

        if(timer && bar_num[i] > 0){for(i in 1:bar_num[i]){pracma::fprintf('-')}}
      }
      if(timer){pracma::fprintf('\n')}
    }
    boot_variance_id = floor(boot_times * c(0.025, 0.975))
    boot_variance_id[boot_variance_id == 0] = 1

    if(get_DE){
      boot_DE_mat = my_sort_mat(boot_DE_mat)
      result$DE$boot_lower = boot_DE_mat[boot_variance_id[1], ]
      result$DE$boot_upper = boot_DE_mat[boot_variance_id[2], ]
    }
    if(get_IE){
      boot_IE_mat = my_sort_mat(boot_IE_mat)
      result$IE$boot_lower = boot_IE_mat[boot_variance_id[1], ]
      result$IE$boot_upper = boot_IE_mat[boot_variance_id[2], ]
    }
  }
  if(plot_result){tryCatch(plot_CHH2020(result), error = function(msg){
    print('Something wrong with the plot function. Please tell me.')
    return(NULL)
  })
  }
  if(sen_ana){tryCatch(plot_sen_ana(result), error = function(msg){
    print('Something wrong with the plot function. Please tell me.')
    return(NULL)
  })
  }
  set.seed(Sys.time())
  return(result)
}

#' Semiparametric Causal Mediation Modeling of Semi-Competing Risks
#'
#' A function that simulates different situations.
#'
#' @param simulation_type Choices are \code{1} (unbiasedness) or \code{2} (coverage).
#' @param hypo Choices are \code{null} or \code{alter}.
#' @param sample_size a number only effective when \code{simulation_type = 2}. Default is \code{1000}.
#' @param repeat_size Default is \code{1000}.
#' @param num_of_cores the number of cores assigned. Default is \code{1}.
#' @param timer will show the progress. Default is \code{TRUE}.
#' @return \code{simulation} returns 6 plots if \code{simulation_type == 1}; a data frame containing coverage if \code{simulation_type == 2}.
#' @export
simulation = function(simulation_type, hypo, sample_size = 1000, repeat_size = 1000, num_of_cores = 1, timer = TRUE, get_variance = c('a', 'b')){
  if(simulation_type == 1){
    plot_successful = unbiasedness(hypo = hypo, sample_size = 1000, repeat_size = repeat_size, num_of_cores = num_of_cores, timer = timer)
    return(NULL)
  }else if(simulation_type == 2){
    coverage_rate = coverage(hypo = hypo, sample_size = sample_size, repeat_size = repeat_size, num_of_cores = num_of_cores, timer = timer, get_variance = get_variance)
    return(coverage_rate)
  }else{
    warning("Unreconginized simulation type.")
    return(NULL)
  }
}
unbiasedness = function(hypo, sample_size, repeat_size, num_of_cores = 1, timer = TRUE){
  result_DE = list()
  result_DE$FF = vector(mode = "list", length = repeat_size)
  result_DE$TF = vector(mode = "list", length = repeat_size)
  result_DE$TT = vector(mode = "list", length = repeat_size)

  result_IE = list()
  result_IE$FF = vector(mode = "list", length = repeat_size)
  result_IE$TF = vector(mode = "list", length = repeat_size)
  result_IE$TT = vector(mode = "list", length = repeat_size)

  min_max_DE = c(0, 0)
  min_max_IE = c(0, 0)

  if(num_of_cores > 1){
    require(foreach)

    ## num_of_cores set-up
    cores = num_of_cores
    cl = snow::makeCluster(cores[1])
    my_functions = c("CASCR", "generate_df", "generate_df2", "alternative_z_1_2", "downsample_func", "data_preprocess", "df_shift_to_cal_level", "do_sen_ana", "estimate_alpha", "estimate_effect", "form_matrix", "get_alpha_variance", "get_beta_variance", "get_counterfactual_hazard", "get_pd", "get_position", "compute_variance", "inv_coxinformation", "make_small", "my_basehaz", "my_eva_fun", "my_sort_mat", "mycoxph", "my_rep_row")
    snow::clusterExport(cl, my_functions)
    doSNOW::registerDoSNOW(cl)
    pb = txtProgressBar(max = repeat_size, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)

    i = 1
    result_now = foreach(i = 1:repeat_size, .options.snow = opts, .combine = 'c', .export = my_functions)%dopar%{
      df_FF = generate_df(sample_size, repeat_size = 1, hypo, confounder = F, calibration = F, myseed = i)
      df_TF = generate_df(sample_size, repeat_size = 1, hypo, confounder = T, calibration = F, myseed = i)
      df_TT = generate_df(sample_size, repeat_size = 1, hypo, confounder = T, calibration = T, myseed = i)

      result_FF = CASCR(df_FF[[1]], get_variance = NULL, timer = FALSE, intervention = c(2, 1))
      result_TF = CASCR(df_TF[[1]], get_variance = NULL, timer = FALSE, intervention = c(2, 1))
      result_TT = CASCR(df_TT[[1]], get_variance = NULL, timer = FALSE, intervention = c(2, 1))
      result_now = list(result_FF = result_FF, result_TF = result_TF, result_TT = result_TT)
      gc()
      return(list(result_now))
    }
    snow::stopCluster(cl)
    pracma::fprintf('\n')

    for(i in 1:repeat_size){
      result_DE$FF[[i]]$time   = result_now[[i]][[1]]$DE$time
      result_DE$FF[[i]]$effect = result_now[[i]][[1]]$DE$effect
      result_DE$FF[[i]]$sick   = result_now[[i]][[1]]$cox_b1$cum_haz$time
      result_DE$TF[[i]]$time   = result_now[[i]][[2]]$DE$time
      result_DE$TF[[i]]$effect = result_now[[i]][[2]]$DE$effect
      result_DE$TF[[i]]$sick   = result_now[[i]][[2]]$cox_b1$cum_haz$time
      result_DE$TT[[i]]$time   = result_now[[i]][[3]]$DE$time
      result_DE$TT[[i]]$effect = result_now[[i]][[3]]$DE$effect
      result_DE$TT[[i]]$sick   = result_now[[i]][[3]]$cox_b1$cum_haz$time

      result_IE$FF[[i]]$time   = result_now[[i]][[1]]$IE$time
      result_IE$FF[[i]]$effect = result_now[[i]][[1]]$IE$effect
      result_IE$FF[[i]]$sick   = result_now[[i]][[1]]$cox_b1$cum_haz$time
      result_IE$TF[[i]]$time   = result_now[[i]][[2]]$IE$time
      result_IE$TF[[i]]$effect = result_now[[i]][[2]]$IE$effect
      result_IE$TF[[i]]$sick   = result_now[[i]][[2]]$cox_b1$cum_haz$time
      result_IE$TT[[i]]$time   = result_now[[i]][[3]]$IE$time
      result_IE$TT[[i]]$effect = result_now[[i]][[3]]$IE$effect
      result_IE$TT[[i]]$sick   = result_now[[i]][[3]]$cox_b1$cum_haz$time

      min_max_DE[1] = min(min_max_DE[1], result_DE$FF[[i]]$effect, result_DE$TF[[i]]$effect, result_DE$TT[[i]]$effect)
      min_max_DE[2] = max(min_max_DE[2], result_DE$FF[[i]]$effect, result_DE$TF[[i]]$effect, result_DE$TT[[i]]$effect)

      min_max_IE[1] = min(min_max_IE[1], result_IE$FF[[i]]$effect, result_IE$TF[[i]]$effect, result_IE$TT[[i]]$effect)
      min_max_IE[2] = max(min_max_IE[2], result_IE$FF[[i]]$effect, result_IE$TF[[i]]$effect, result_IE$TT[[i]]$effect)
    }
  }else{
    ## fetch basic parameter
    if(timer){
      space = 100
      pracma::fprintf('|         10       20        30        40        50        60        70        80        90    100 |\n')
      loop_count = 1:repeat_size
      counter_total = repeat_size
      cum_bar_num = my_eva_fun(list(1:space, 1:space / space * counter_total), loop_count, rule = '0')
      bar_num = diff(c(0, cum_bar_num))
    }
    i = 1
    for(i in 1:repeat_size){
      # print(i)
      df_FF = generate_df(sample_size, repeat_size = 1, hypo, confounder = F, calibration = F, myseed = i)
      df_TF = generate_df(sample_size, repeat_size = 1, hypo, confounder = T, calibration = F, myseed = i)
      df_TT = generate_df(sample_size, repeat_size = 1, hypo, confounder = T, calibration = T, myseed = i)

      result_FF = CASCR(df_FF[[1]], get_variance = NULL, timer = FALSE, intervention = c(2, 1))
      result_TF = CASCR(df_TF[[1]], get_variance = NULL, timer = FALSE, intervention = c(2, 1))
      result_TT = CASCR(df_TT[[1]], get_variance = NULL, timer = FALSE, intervention = c(2, 1))

      result_DE$FF[[i]]$time   = result_FF$DE$time
      result_DE$FF[[i]]$effect = result_FF$DE$effect
      result_DE$FF[[i]]$sick   = result_FF$cox_b1$cum_haz$time
      result_DE$TF[[i]]$time   = result_TF$DE$time
      result_DE$TF[[i]]$effect = result_TF$DE$effect
      result_DE$TF[[i]]$sick   = result_TF$cox_b1$cum_haz$time
      result_DE$TT[[i]]$time   = result_TT$DE$time
      result_DE$TT[[i]]$effect = result_TT$DE$effect
      result_DE$TT[[i]]$sick   = result_TT$cox_b1$cum_haz$time

      result_IE$FF[[i]]$time   = result_FF$IE$time
      result_IE$FF[[i]]$effect = result_FF$IE$effect
      result_IE$FF[[i]]$sick   = result_FF$cox_b1$cum_haz$time
      result_IE$TF[[i]]$time   = result_TF$IE$time
      result_IE$TF[[i]]$effect = result_TF$IE$effect
      result_IE$TF[[i]]$sick   = result_TF$cox_b1$cum_haz$time
      result_IE$TT[[i]]$time   = result_TT$IE$time
      result_IE$TT[[i]]$effect = result_TT$IE$effect
      result_IE$TT[[i]]$sick   = result_TT$cox_b1$cum_haz$time

      min_max_DE[1] = min(min_max_DE[1], result_DE$FF[[i]]$effect, result_DE$TF[[i]]$effect, result_DE$TT[[i]]$effect)
      min_max_DE[2] = max(min_max_DE[2], result_DE$FF[[i]]$effect, result_DE$TF[[i]]$effect, result_DE$TT[[i]]$effect)

      min_max_IE[1] = min(min_max_IE[1], result_IE$FF[[i]]$effect, result_IE$TF[[i]]$effect, result_IE$TT[[i]]$effect)
      min_max_IE[2] = max(min_max_IE[2], result_IE$FF[[i]]$effect, result_IE$TF[[i]]$effect, result_IE$TT[[i]]$effect)
      if(timer && bar_num[i] > 0){for(i in 1:bar_num[i]){pracma::fprintf('-')}}
    }
    if(timer){pracma::fprintf('\n')}
  }

  width = 500
  height = 500
  # png(file = paste("/Users/js/Desktop/CHH2020/bias_DE_", hypo, "_no_conf.png", sep = ''), width = width, height = height)
  true_DE = alternative_z_1_2(hypo, effect = 'DE', confounder = F, intervention = c(2, 1), time_by = 1e-2)
  plot_successful = plot_unbiasedness(result_DE$FF, true_DE, ylim = min_max_DE, hypo, effect = 'DE', confounder = F, calibration = F)
  # dev.off()
  # png(file = paste("/Users/js/Desktop/CHH2020/bias_DE_", hypo, "_unadj_conf.png", sep = ''), width = width, height = height)
  true_DE = alternative_z_1_2(hypo, effect = 'DE', confounder = T, intervention = c(2, 1), time_by = 1e-2)
  plot_successful = plot_unbiasedness(result_DE$TF, true_DE, ylim = min_max_DE, hypo, effect = 'DE', confounder = T, calibration = F)
  # dev.off()
  # png(file = paste("/Users/js/Desktop/CHH2020/bias_DE_", hypo, "_adj_conf.png", sep = ''), width = width, height = height)
  true_DE = alternative_z_1_2(hypo, effect = 'DE', confounder = T, intervention = c(2, 1), time_by = 1e-2)
  plot_successful = plot_unbiasedness(result_DE$TT, true_DE, ylim = min_max_DE, hypo, effect = 'DE', confounder = T, calibration = T)
  # dev.off()

  if(hypo == "null"){
    min_max_IE = c(-0.5, 0.5)
  }else{
    min_max_IE = c(-0.6, 0.3)
  }

  true_IE = alternative_z_1_2(hypo, effect = 'IE', confounder = F, intervention = c(2, 1), time_by = 1e-2)
  # png(file = paste("/Users/js/Desktop/CHH2020/bias_IE_", hypo, "_no_conf.png", sep = ''), width = width, height = height)
  plot_successful = plot_unbiasedness(result_IE$FF, true_IE, ylim = min_max_IE, hypo, effect = 'IE', confounder = F, calibration = F)
  # dev.off()
  true_IE = alternative_z_1_2(hypo, effect = 'IE', confounder = T, intervention = c(2, 1), time_by = 1e-2)
  # png(file = paste("/Users/js/Desktop/CHH2020/bias_IE_", hypo, "_unadj_conf.png", sep = ''), width = width, height = height)
  plot_successful = plot_unbiasedness(result_IE$TF, true_IE, ylim = min_max_IE, hypo, effect = 'IE', confounder = T, calibration = F)
  # dev.off()
  true_IE = alternative_z_1_2(hypo, effect = 'IE', confounder = T, intervention = c(2, 1), time_by = 1e-2)
  # png(file = paste("/Users/js/Desktop/CHH2020/bias_IE_", hypo, "_adj_conf.png", sep = ''), width = width, height = height)
  plot_successful = plot_unbiasedness(result_IE$TT, true_IE, ylim = min_max_IE, hypo, effect = 'IE', confounder = T, calibration = T)
  # dev.off()
  return(TRUE)
}
coverage = function(hypo, sample_size, repeat_size, num_of_cores, timer = TRUE, get_variance = c('a', 'b')){
  if(hypo == 'alter'){
    true_DE = alternative_z_1_2(hypo, effect = 'DE', confounder = F, intervention = c(2, 1))
    true_IE = alternative_z_1_2(hypo, effect = 'IE', confounder = F, intervention = c(2, 1))
  }
  if(num_of_cores > 1){
    require(foreach)

    ## num_of_cores set-up
    cores = num_of_cores
    cl = snow::makeCluster(cores[1])
    my_functions = c("CASCR", "generate_df", "generate_df2", "alternative_z_1_2", "downsample_func", "data_preprocess", "df_shift_to_cal_level", "do_sen_ana", "estimate_alpha", "estimate_effect", "form_matrix", "get_alpha_variance", "get_beta_variance", "get_counterfactual_hazard", "get_pd", "get_position", "compute_variance", "inv_coxinformation", "make_small", "my_basehaz", "my_eva_fun", "my_sort_mat", "mycoxph", "my_rep_row")
    snow::clusterExport(cl, my_functions)
    doSNOW::registerDoSNOW(cl)
    pb = txtProgressBar(max = repeat_size, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)

    i = 1; myunit = 'raw'; variance_method = 'new'
    result_now = foreach(i = 1:repeat_size, .options.snow = opts, .combine = 'c', .export = my_functions)%dopar%{
      df1 = generate_df(sample_size, repeat_size = 1, hypo, confounder = F, calibration = F, myseed = i)[[1]]
      result_1 = CASCR(df1, get_variance = get_variance, timer = FALSE, intervention = c(2, 1), myunit = myunit, variance_method = variance_method)
      if(is.null(result_1$DE$boot_lower)){
        result_1$DE$boot_lower = result_1$DE$asym_lower
        result_1$DE$boot_upper = result_1$DE$asym_upper
        result_1$IE$boot_lower = result_1$IE$asym_lower
        result_1$IE$boot_upper = result_1$IE$asym_upper
      }
      ind = floor(length(result_1$DE$time) * c(0.2, 0.4, 0.5, 0.6, 0.8))
      result_1 = list(time = result_1$DE$time[ind],
                      DE = data.frame(asym_lower = result_1$DE$asym_lower[ind], asym_upper = result_1$DE$asym_upper[ind], boot_lower = result_1$DE$boot_lower[ind], boot_upper = result_1$DE$boot_upper[ind]),
                      IE = data.frame(asym_lower = result_1$IE$asym_lower[ind], asym_upper = result_1$IE$asym_upper[ind], boot_lower = result_1$IE$boot_lower[ind], boot_upper = result_1$IE$boot_upper[ind]))

      if(hypo == 'null'){
        df2 = generate_df2(sample_size, myseed = i)
        result_2 = CASCR(df2, get_variance = get_variance, timer = FALSE, intervention = c(2, 1), myunit = myunit, variance_method = variance_method)
        if(is.null(result_2$DE$boot_lower)){
          result_2$DE$boot_lower = result_2$DE$asym_lower
          result_2$DE$boot_upper = result_2$DE$asym_upper
          result_2$IE$boot_lower = result_2$IE$asym_lower
          result_2$IE$boot_upper = result_2$IE$asym_upper
        }
        ind = floor(length(result_2$IE$time) * c(0.2, 0.4, 0.5, 0.6, 0.8))
        result_2 = list(DE = data.frame(asym_lower = result_2$DE$asym_lower[ind], asym_upper = result_2$DE$asym_upper[ind], boot_lower = result_2$DE$boot_lower[ind], boot_upper = result_2$DE$boot_upper[ind]),
                        IE = data.frame(asym_lower = result_2$IE$asym_lower[ind], asym_upper = result_2$IE$asym_upper[ind], boot_lower = result_2$IE$boot_lower[ind], boot_upper = result_2$IE$boot_upper[ind]))
        result_now = list(DE = result_1$DE, IE = result_1$IE, IE2 = result_2$IE)
      }else{
        true_DE_now = approx(true_DE$time, true_DE$hazard, xout = result_1$time, rule = 2)$y
        true_IE_now = approx(true_IE$time, true_IE$hazard, xout = result_1$time, rule = 2)$y
        result_1$DE = result_1$DE - true_DE_now
        result_1$IE = result_1$IE - true_IE_now
        result_now = list(DE = result_1$DE, IE = result_1$IE)
      }
      gc()
      return(list(result_now))
    }
    snow::stopCluster(cl)
    pracma::fprintf('\n')

  }else{
    result_now = vector(mode = 'list', length = repeat_size)
    ## fetch basic parameter
    if(timer){
      space = 100
      pracma::fprintf('|         10       20        30        40        50        60        70        80        90    100 |\n')
      loop_count = 1:repeat_size
      counter_total = repeat_size
      cum_bar_num = my_eva_fun(list(1:space, 1:space / space * counter_total), loop_count, rule = '0')
      bar_num = diff(c(0, cum_bar_num))
    }
    i = 1; myunit = 'raw'; variance_method = 'new'
    for(i in 1:repeat_size){
      df1 = generate_df(sample_size, repeat_size = 1, hypo, confounder = F, calibration = F, myseed = i)[[1]]
      result_1 = CASCR(df1, get_variance = get_variance, timer = FALSE, intervention = c(2, 1), myunit = myunit, variance_method = variance_method)
      if(is.null(result_1$DE$boot_lower)){
        result_1$DE$boot_lower = result_1$DE$asym_lower
        result_1$DE$boot_upper = result_1$DE$asym_upper
        result_1$IE$boot_lower = result_1$IE$asym_lower
        result_1$IE$boot_upper = result_1$IE$asym_upper
      }
      ind = floor(length(result_1$DE$time) * c(0.2, 0.4, 0.5, 0.6, 0.8))
      result_1 = list(time = result_1$DE$time[ind],
                      DE = data.frame(asym_lower = result_1$DE$asym_lower[ind], asym_upper = result_1$DE$asym_upper[ind], boot_lower = result_1$DE$boot_lower[ind], boot_upper = result_1$DE$boot_upper[ind]),
                      IE = data.frame(asym_lower = result_1$IE$asym_lower[ind], asym_upper = result_1$IE$asym_upper[ind], boot_lower = result_1$IE$boot_lower[ind], boot_upper = result_1$IE$boot_upper[ind]))

      if(hypo == 'null'){
        df2 = generate_df2(sample_size, myseed = i)
        result_2 = CASCR(df2, get_variance = get_variance, timer = FALSE, intervention = c(2, 1), myunit = myunit, variance_method = variance_method)
        if(is.null(result_2$DE$boot_lower)){
          result_2$DE$boot_lower = result_2$DE$asym_lower
          result_2$DE$boot_upper = result_2$DE$asym_upper
          result_2$IE$boot_lower = result_2$IE$asym_lower
          result_2$IE$boot_upper = result_2$IE$asym_upper
        }
        ind = floor(length(result_2$IE$time) * c(0.2, 0.4, 0.5, 0.6, 0.8))
        result_2 = list(DE = data.frame(asym_lower = result_2$DE$asym_lower[ind], asym_upper = result_2$DE$asym_upper[ind], boot_lower = result_2$DE$boot_lower[ind], boot_upper = result_2$DE$boot_upper[ind]),
                        IE = data.frame(asym_lower = result_2$IE$asym_lower[ind], asym_upper = result_2$IE$asym_upper[ind], boot_lower = result_2$IE$boot_lower[ind], boot_upper = result_2$IE$boot_upper[ind]))
        result_now[[i]] = list(DE = result_1$DE, IE = result_1$IE, IE2 = result_2$IE)
      }else{
        true_DE_now = approx(true_DE$time, true_DE$hazard, xout = result_1$time, rule = 2)$y
        true_IE_now = approx(true_IE$time, true_IE$hazard, xout = result_1$time, rule = 2)$y
        result_1$DE = result_1$DE - true_DE_now
        result_1$IE = result_1$IE - true_IE_now
        result_now[[i]] = list(DE = result_1$DE, IE = result_1$IE)
      }
      if(timer && bar_num[i] > 0){for(i in 1:bar_num[i]){pracma::fprintf('-')}}
    }
    if(timer){pracma::fprintf('\n')}
  }
  BootVariance = sum(c('b', 'B', 'boot', 'bootstrap', 'bootstrapping', 'Boot', 'Bootstrap', 'Bootstrapping') %in% get_variance) > 0

  if(BootVariance){
    DE_coverage = data.frame(asym = rep(0, 5), boot = rep(0, 5))
    IE_coverage = data.frame(asym = rep(0, 5), boot = rep(0, 5))
    if(hypo == "null"){IE2_coverage = data.frame(asym = rep(0, 5), boot = rep(0, 5))}
  }else{
    DE_coverage = data.frame(asym = rep(0, 5))
    IE_coverage = data.frame(asym = rep(0, 5))
    if(hypo == "null"){IE2_coverage = data.frame(asym = rep(0, 5))}
  }

  for(i in 1:repeat_size){
    DE_coverage$asym = DE_coverage$asym + ((result_now[[i]]$DE$asym_lower * result_now[[i]]$DE$asym_upper) < 0)
    if(BootVariance){DE_coverage$boot = DE_coverage$boot + ((result_now[[i]]$DE$boot_lower * result_now[[i]]$DE$boot_upper) < 0)}

    IE_coverage$asym = IE_coverage$asym + ((result_now[[i]]$IE$asym_lower * result_now[[i]]$IE$asym_upper) < 0)
    if(BootVariance){IE_coverage$boot = IE_coverage$boot + ((result_now[[i]]$IE$boot_lower * result_now[[i]]$IE$boot_upper) < 0)}

    if(hypo == "null"){
      IE2_coverage$asym = IE2_coverage$asym + ((result_now[[i]]$IE2$asym_lower * result_now[[i]]$IE2$asym_upper) < 0)
      if(BootVariance){IE2_coverage$boot = IE2_coverage$boot + ((result_now[[i]]$IE2$boot_lower * result_now[[i]]$IE2$boot_upper) < 0)}
    }
  }
  if(hypo == "null"){
    return(list(DE = DE_coverage/repeat_size, IE1 = IE_coverage/repeat_size, IE2 = IE2_coverage/repeat_size))
  }else{
    return(list(DE = DE_coverage/repeat_size, IE = IE_coverage/repeat_size))
  }
}

## simulation
#' @export
generate_df = function(sample_size, repeat_size, hypo, confounder, calibration, myseed = 1){
  df_all = vector(mode = 'list', length = repeat_size)
  set.seed(myseed)
  alpha1Z = 0.25 * (hypo == 'alter')
  alpha2Z = 0.25 * (hypo == 'alter')
  intersection = 0
  calibration = calibration & confounder

  alphaX = 1 * confounder
  X = c(rep(0, sample_size/2), rep(1, sample_size/2)) * confounder
  Z = ((X + rnorm(sample_size)) > 0.5) + 1
  for(counter in 1:repeat_size){
    set.seed(counter + myseed)

    T1 = rweibull(sample_size, scale = exp(intersection + alpha1Z * Z + alphaX * X), shape = 1)
    T2 = T1 + 0.5 * rweibull(sample_size, scale = exp(alpha2Z * Z + alphaX * X), shape = 1)
    C = rweibull(sample_size, scale = 2, shape = 5)

    d2 = T2 < C
    T2 = pmin(T2, C)
    d1 = T1 < T2
    T1 = pmin(T1, T2)

    if(calibration){df_all[[counter]] = data.frame(T1 = T1, T2 = T2, d1 = d1, d2 = d2, Z = Z, X = X)}
    if(!calibration){df_all[[counter]] = data.frame(T1 = T1, T2 = T2, d1 = d1, d2 = d2, Z = Z)}
  }
  return(df_all)
}
#' @export
generate_df2 = function(sample_size, myseed = 1){
  set.seed(myseed)
  Z = c(rep(1, sample_size/2), rep(2, sample_size/2))
  T1 = runif(sample_size, 0, 2)
  picked_index = runif(sample_size) < 0.75
  T1[Z == 2 & picked_index] = runif(sum(Z == 2 & picked_index), 1.5, 2)

  # T1 = rbeta(sample_size, shape1 = 1, shape2 = 2)
  # T1[Z == 2] = rbeta(sum(Z == 2), shape1 = 2, shape2 = 1)

  T2 = 2 * rbeta(sample_size, shape1 = 2, shape2 = 1)
  C = rweibull(sample_size, scale = 2, shape = 5)

  d2 = T2 < C
  T2 = pmin(T2, C)
  d1 = T1 < T2
  T1 = pmin(T1, T2)

  df = data.frame(T1 = T1, T2 = T2, d1 = d1, d2 = d2, Z = Z)
  return(df)
}
#' @export
alternative_z_1_2 = function(hypo, effect, confounder, intervention, time_by = 5e-4){
  #### you can change the parameters
  alpha1Z = 0.25
  alpha2Z = 0.25
  alphaX = 1
  intersection = 0

  #### don't touch me
  tstart = 0
  tend = 4
  t = seq(tstart, tend, by = time_by)
  diff_t = t[2] - t[1]

  alpha1Z = alpha1Z * (hypo == 'alter')
  alpha2Z = alpha2Z * (hypo == 'alter')
  alphaX = alphaX * confounder * 0.5

  #### (2, 2)
  z_a = intervention[1]
  z_b = ifelse(effect == 'DE', intervention[2], intervention[1])
  a1zb = exp(intersection + alpha1Z * z_b + alphaX)
  ca2zb = 0.5 * exp(alpha2Z * z_b + alphaX)
  ca2za = 0.5 * exp(alpha2Z * z_a + alphaX)
  w0zb = (1 - sdprisk::phypoexp(t, 1/a1zb)) / (1 - sdprisk::phypoexp(t, c(1/a1zb, 1/ca2zb * (1 + 1e-7))))
  w1zb = 1 - w0zb
  case1 = cumsum(w1zb) * diff_t / ca2za

  #### (2, 1)
  z_a = ifelse(effect == 'DE', intervention[2], intervention[1])
  z_b = intervention[2]
  a1zb = exp(intersection + alpha1Z * z_b + alphaX)
  ca2zb = 0.5 * exp(alpha2Z * z_b + alphaX)
  ca2za = 0.5 * exp(alpha2Z * z_a + alphaX)
  w0zb = (1 - sdprisk::phypoexp(t, 1/a1zb)) / (1 - sdprisk::phypoexp(t, c(1/a1zb, 1/ca2zb * (1 + 1e-7))))
  w1zb = 1 - w0zb
  case2 = cumsum(w1zb) * diff_t / ca2za

  #### (2, 2) - (2, 1)
  cumhaz = data.frame(hazard = case1 - case2, time = t)

  return(cumhaz)
}

## data related
#' @export
downsample_func = function(vec, downsample){
  if(length(vec) %% downsample == 1){
    return(vec[seq(1, length(vec), by = downsample)])
  }else{
    return(vec[c(seq(1, length(vec), by = downsample), length(vec))])
  }
}
#' @export
data_preprocess = function(df, myunit, downsample){
  set.seed(1)
  ## remove unreasonable columns
  if(sum(is.na(df)) > 0){stop('NA exists.')}
  col_var = apply(as.matrix(df[5:dim(df)[2]]), 2, var) < .Machine$double.eps
  if(col_var[1]){stop("All exposures are identical.")}
  if(sum(col_var) > 0){
    col_var_logi = which(col_var)
    msg = paste0('column(s) ', colnames(df)[col_var_logi + 4], ' have variance = 0, which will be removed.')
    warning(msg, immediate. = TRUE)
    df = data.frame(df[, 1:4], df[, 5 + col_var_logi])
  }

  ## rename columns and adjust data
  colnames(df)[1:4] = c('T1', 'T2', 'd1', 'd2')
  df$d1 = as.logical(df$d1)
  df$d2 = as.logical(df$d2)

  df = df[!(df$T1 == 0 & df$T2 == 0 & df$d1 == 0 & df$d2 == 0), ]
  df$T1[df$T1 == 0 & df$T2 == 0 & df$d1 == 1] = 1
  df$T1[df$T1 == 0 & df$T2 == 0 & df$d1 == 0 & df$d2 == 1] = 2
  df$T2[df$T1 == 0 & df$T2 == 0 & df$d2 == 1] = 2
  negative_shift = abs(apply(df[df$T1 < 0 | df$T2 < 0, c(1, 2)], 1, min))
  df$T1[df$T1 < 0 | df$T2 < 0] = df$T1[df$T1 < 0 | df$T2 < 0] + negative_shift
  df$T2[df$T1 < 0 | df$T2 < 0] = df$T2[df$T1 < 0 | df$T2 < 0] + negative_shift

  # df$d1[df$d1 & (df$T1 == df$T2)] = FALSE
  df$T2[df$d1 & (df$T1 == df$T2)] = df$T2[df$d1 & (df$T1 == df$T2)] + 1

  # df$d1[df$d1 == 0 & df$T1 < df$T2] = TRUE
  df$T1[df$d1 == 0 & df$T1 < df$T2] = df$T2[df$d1 == 0 & df$T1 < df$T2]

  ## reunit
  if(myunit == 'raw'){
    T1 = df$T1
    T2 = df$T2
    if(downsample > 1){
      if(abs(downsample - round(downsample)) > .Machine$double.eps){
        downsample = ceiling(downsample)
        warning(paste0('\'downsample\' should be an integer. It will be rounded to ', downsample, '.'), immediate. = TRUE)
      }
      ## get real time for two Coxs
      sort_T2_obs = sort(T2[df$d2 == 1], index.return = TRUE)
      sort_D1 = (df$d1[df$d2 == 1])[sort_T2_obs$ix]

      downsampled_T2 = downsample_func(sort_T2_obs$x, downsample)
      downsampled_D11 = diff(c(0, approx(x = sort_T2_obs$x, y = cumsum(sort_D1), xout = downsampled_T2, method = 'constant', ties = 'max')$y)) > 0
      downsampled_D10 = diff(c(0, approx(x = sort_T2_obs$x, y = cumsum(!sort_D1), xout = downsampled_T2, method = 'constant', ties = 'max')$y)) > 0

      unique_T2 = unique(downsampled_T2)
      b0_time = unique(downsampled_T2[downsampled_D10 == 1])
      b1_time = unique(downsampled_T2[downsampled_D11 == 1])
    }
  }else{
    T1 = floor(df$T1/myunit) * myunit
    T1[T1 == 0] = myunit
    T2 = floor(df$T2/myunit) * myunit
    T2[T2 == 0] = myunit
    if(downsample != 1){
      downsample = 1
      warning('\'downsample\' will be ignored when \'myunit\' is on.', immediate. = TRUE)
    }
  }
  if(downsample == 1){
    ## get real time for two Coxs
    unique_T2 = unique(sort(T2[df$d2 == 1]))
    b0_time = unique(sort(T1[df$d2 & (df$d1 == 0)]))
    b1_time = unique(sort(T2[df$d2 & (df$d1 == 1)]))
  }

  ## randomly and slightly move time points
  m = dim(df)[1]
  t_diff = diff(sort(df$T2))
  t_diff = min(t_diff[t_diff > 0])
  perturbation = t_diff/100 * (runif(m) + runif(m))
  df$T1 = df$T1 + perturbation
  df$T2 = df$T2 + perturbation

  ## rearrange df such that [df$primary_observed = 1; df$primary_observed = 0]
  primary_observed_num = sum(df$d2)
  df = rbind(df[df$d2==1, ], df[df$d2==0, ])

  ## sort both two part by death time and censoring time respectively.
  # Not censored part is sorted by death time.
  dataNC = df[1:primary_observed_num, ]
  dataNC = dataNC[order(dataNC$T2), ]

  # Censored part (if exists) is sorted by censoring time.
  if(primary_observed_num + 1 <= m){
    dataC = df[(primary_observed_num+1):m, ]
    dataC = dataC[order(dataC$T2), ]
    df = rbind(dataNC, dataC)
  }else{
    df = dataNC
  }
  df = list(df = df, unique_T2 = unique_T2, b0_time = b0_time, b1_time = b1_time)
  return(df)
}
#' @export
df_shift_to_cal_level = function(df, cal_level){
  num_cal = dim(df)[2] - 5
  if(num_cal == 0){
    cal_level = NULL
  }else if(is.character(cal_level)){
    cal_level = sapply(df, cal_level)[6:(5 + num_cal)]
  }else if(length(cal_level) != num_cal){
    warning('length of cal_level is wrong. Default (median) is used.', immediate. = TRUE)
    cal_level = sapply(df, median)[6:(5 + num_cal)]
  }
  if(num_cal > 0){
    df[, 6:(5 + num_cal)] = df[, 6:(5 + num_cal)] - my_rep_row(cal_level, dim(df)[1])
    cal_level = rep(0, dim(df)[2] - 5)
  }
  return(list(df = df, cal_level = cal_level))
}

## others
#' @export
make_small = function(cox_b, unique_T2){
  # cox_b = cox_b0
  if(length(cox_b) == 1 || length(cox_b) == 0){
    return(cox_b)
  }else{
    all_time_diff = diff(sort(c(cox_b$cum_haz$cum_haz, unique_T2)))
    time_diff = min(all_time_diff[all_time_diff!=0])/3
    new_cum_haz = approx(x = cox_b$cum_haz$time, y = cox_b$cum_haz$cum_haz, xout = unique_T2 + time_diff, rule = 2, method = 'constant')
    cox_b$cum_haz = data.frame(cum_haz = new_cum_haz$y, time = unique_T2)
    return(cox_b)
  }
}
#' @export
my_basehaz = function(time, observed, covariates, cox){
  obs_time = time[observed, 2]

  time_data = cbind(time, observed)
  ## lower-left ## 1 for start; 2 for end;
  rank_time_2 = rank(time_data[time_data[, 3] == 1, 2], ties.method = 'min')
  sort_time_2 = sort(time_data[, 2], index.return = TRUE)
  sorted_time_df_2 = time_data[sort_time_2$ix, ]

  sorted_covariates = as.matrix(covariates[sort_time_2$ix, ])
  exp_tmp = exp(sorted_covariates %*% cox$coefficients)
  tmp_lower_left = (sum(exp_tmp) - c(0, cumsum(exp_tmp)))[-(dim(exp_tmp)[1] + 1)]
  tmp_lower_left = tmp_lower_left[sorted_time_df_2[, 3] == 1]
  tmp2 = tmp_lower_left[rank_time_2]

  sort_time_1 = sort(time_data[, 1], index.return = TRUE)
  sorted_time_df_1 = time_data[sort_time_1$ix, ]
  important_index = get_position(time_data[time_data[, 3] == 1, 2], sort_time_1$x)

  sorted_covariates = as.matrix(covariates[sort_time_1$ix, ])
  exp_tmp = exp(sorted_covariates %*% cox$coefficients)
  tmp_lower_left = (sum(exp_tmp) - c(0, cumsum(exp_tmp))) # use important index
  tmp1 = tmp_lower_left[important_index]

  cum_haz = data.frame(cum_haz = cumsum(1/(tmp2 - tmp1)), time = obs_time)
  return(cum_haz)
}
#' @export
my_rep_row = function(x, n){
  return(matrix(rep(x, each = n), nrow = n))
}
#' @export
form_matrix = function(sd_time_indep, sd_time_mix, sd_time_dep){
  s_dim = dim(sd_time_indep)[1]
  l_dim = s_dim + length(sd_time_dep)
  m = matrix(0, l_dim, l_dim)
  m[1:s_dim, 1:s_dim] = sd_time_indep
  m[1:s_dim, (s_dim+1):l_dim] = t(sd_time_mix)
  m[(s_dim+1):l_dim, 1:s_dim] = sd_time_mix
  m[(s_dim+1):l_dim, (s_dim+1):l_dim] = diag(sd_time_dep)
  return(m)
}
#' @export
my_eva_fun = function(fun1, points, rule = '0', method = 'constant'){
  if(rule == 'no0'){
    return(approx(fun1[[2]], fun1[[1]], points, method = method, rule = 2, ties = max)$y)
  }else if(rule == '0'){
    return(approx(c(0, fun1[[2]]), c(0, fun1[[1]]), points, method = method, rule = 2, ties = max)$y)
  }
}
#' @export
get_position = function(x, y){
  # get_position returns vector z.
  # z_i := min_j(x_i<=y_j)
  # z_i := length(y) + 1, if x_i>y_j for all j
  # y is a ordered sequence. x may or may not be ordered. x and y can have different lengths.
  if(is.unsorted(y)){
    stop('The second argument must be sorted.')
  }
  m = length(x)
  z = rep(0, m)
  now_index = 1
  reach_end = 0
  tmp_x = sort(x)
  position_x = order(x)
  for(i in 1:m){
    if(reach_end == 0){
      while(tmp_x[i] > y[now_index]){
        now_index = now_index + 1
        if(now_index == length(y) + 1){
          reach_end = 1
          break
        }
      }
      z[position_x[i]] = now_index
    }
    if(reach_end == 1){
      z[position_x[i]] = now_index
    }
  }
  return(z)
}

## bootstrap variance
#' @export
my_sort_mat = function(mat){
  mat_na = is.na(mat)
  if(sum(mat_na) == 0){
    mat = apply(mat, 2, sort)
  }else{
    mat = mat[rowSums(mat_na) == 0, ]
    mat = apply(mat, 2, sort)
  }
  return(mat)
}

## estimation
#' @export
estimate_alpha = function(df, cal_level, cox_b1, cox_whole, unique_T2, get_variance, timer, num_of_cores, variance_method, threshold){
  # cox_b1 = small_cox_b1; cox_whole = small_cox_whole

  ## fetch basic parameters
  m = dim(df)[1]
  n_col = dim(df)[2]
  n_covariates = n_col - 4 + 1
  observed_t2 = length(unique_T2)
  AsymVariance = sum(c('a', 'A', 'asym', 'asymptotic', 'asymptotical', 'Asym', 'Asymptotic', 'Asymptotical') %in% get_variance) > 0

  ## preallocation
  # for alpha
  alpha_mat = matrix(0, n_covariates, observed_t2)
  if(AsymVariance){
    n_covariates_cov = n_covariates - length(cal_level)
    alpha_var = matrix(0, n_covariates * observed_t2, n_covariates)
    alpha_cov = matrix(0, n_covariates_cov * observed_t2, n_covariates_cov * observed_t2)
  }

  ## auxiliary
  sort_t2 = sort(df$T2, index.return = TRUE)
  rank_t2 = approx(x = c(0, sort_t2$x), y = 0:m, xout = unique_T2, rule = 2, method = 'constant')$y + 1
  sort_exact_t2 = sort(df$T2[df$d2 == 1])
  rank_exact_t2 = approx(x = c(0, sort_exact_t2), y = 0:length(sort_exact_t2), xout = unique_T2, ties = 'max', rule = 2, method = 'constant')$y + 1
  order_belonged_to = approx(x = rank_exact_t2, y = 1:observed_t2, xout = 1:sum(df$d2), ties = 'max', method = 'constant', rule = 2)$y
  exact_time = list(time = sort_exact_t2, order_belonged_to = order_belonged_to, rank_exact_t2 = rank_exact_t2)

  ## how "large" is the data
  sick_alive = rep(0, observed_t2)
  healthy_alive = rep(0, observed_t2)
  alive = rep(0, observed_t2)

  ## convergence or speed related
  df_alpha_covariates = as.matrix(cbind(rep(1, m), df[sort_t2$ix, 4 + 1:(n_covariates - 1)]))
  df_alpha_time = df[sort_t2$ix, c(1, 3)]
  hard_hreshold = 1e-7
  converged_alpha = rep(FALSE, observed_t2)
  hard_converged_alpha = rep(FALSE, observed_t2)
  t2_index = which(df$d2[sort_t2$ix])
  counter = 0

  if(timer){
    space = 100
    pracma::fprintf('| point estimation 20        30        40        50        60        70        80        90    100 |\n')
    loop_count = 1:observed_t2
    counter_total = observed_t2
    cum_bar_num = my_eva_fun(list(1:space, 1:space / space * counter_total), loop_count, rule = '0')
    bar_num = diff(c(0, cum_bar_num))
  }
  ## estimation for alpha
  for(counter in 1:observed_t2){
    ## counter
    T2 = unique_T2[counter]
    i = rank_t2[counter]
    index_now = ((counter - 1) * n_covariates + 1) : (counter * n_covariates)

    ## split data and get basic information
    sub_Y = (df_alpha_time$T1[i:m] < T2) & (df_alpha_time$d1[i:m] == 1)
    sick_alive[counter] = sum(sub_Y)
    alive[counter] = m - i + 1
    healthy_alive[counter] = alive[counter] - sick_alive[counter]

    if(sum(sub_Y) < 2){
      if(timer && bar_num[counter] > 0){for(k in 1:bar_num[counter]){pracma::fprintf('-')}}
      next
    }
    sub_x = df_alpha_covariates[i:m, ]

    ## estimation for coefficients
    if(counter == 1){
      sub_glm = glm.fit(x = sub_x, y = sub_Y, family = binomial(), intercept = FALSE)
    }else if(!hard_converged_alpha[counter - 1]){
      sub_glm = glm.fit(x = sub_x, y = sub_Y, family = binomial(), intercept = FALSE)
    }else{
      sub_glm = glm.fit(x = sub_x, y = sub_Y, family = binomial(), intercept = FALSE, start = sub_glm$coefficients, etastart = sub_x%*%alpha_vec)
      if(!sub_glm$converged){
        sub_glm = glm.fit(x = sub_x, y = sub_Y, family = binomial(), intercept = FALSE)
      }
    }

    coeff_NA = sum(is.na(sub_glm$coefficients))
    prob_min = min(sub_glm$fitted.values)
    prob_max = max(sub_glm$fitted.values)

    if(prob_min < threshold || prob_max > 1 - threshold || coeff_NA > 0){
      converged_alpha[counter] = FALSE
    }else{
      if(!(prob_min < threshold || prob_max > 1 - threshold)){hard_converged_alpha[counter] = TRUE}
      converged_alpha[counter] = TRUE
      alpha_vec = sub_glm$coefficients
      alpha_mat[, counter] = alpha_vec
    }

    ## estimation of variance
    if(AsymVariance){
      if(converged_alpha[counter]){
        prob = sub_glm$fitted.values
        alpha_var[index_now, ] = solve(crossprod(sub_x * (prob * (1 - prob)), sub_x))
      }
      index_now = index_now + n_covariates
    }
    if(timer && bar_num[counter] > 0){for(k in 1:bar_num[counter]){pracma::fprintf('-')}}
    # print(c(counter, prob_min, prob_max))
  }
  if(timer){pracma::fprintf('\n')}
  # refill alpha who does not converge
  refill_index = approx(x = which(converged_alpha), y = which(converged_alpha), xout = 1:observed_t2, method = 'constant', rule = 2)$y
  alpha_mat = alpha_mat[, refill_index]
  if(AsymVariance){
    for(i in which(!converged_alpha)){
      j = refill_index[i]
      sub_x = df_alpha_covariates[t2_index[i]:m, ]
      prob = 1 / (1 + exp(-as.vector(sub_x %*% alpha_mat[, i])))
      alpha_var[((i - 1) * n_covariates + 1) : (i * n_covariates), ] = tryCatch(solve(crossprod(sub_x * (prob * (1 - prob)), sub_x)), error = function(msg){return(alpha_var[((i - 2) * n_covariates + 1) : ((i - 1) * n_covariates), ])})
    }
  }
  # make some dataframes
  sick_alive = data.frame(time = unique_T2, number = sick_alive)
  alive = data.frame(time = unique_T2, number = alive)
  healthy_alive = data.frame(time = unique_T2, number = healthy_alive)
  converged_alpha = data.frame(time = unique_T2, converged = converged_alpha)

  ## estimation of covariance
  if(AsymVariance){
    alpha_var_tmp = as.matrix(alpha_var[, 1:n_covariates_cov])

    if(variance_method == "new"){
      cum_haz_1_fulltime = approx(x = cox_b1$cum_haz$time, y = cox_b1$cum_haz$cum_haz, xout = unique_T2, method = 'linear', yleft = 0, rule = 2)$y
      survival1_cov = exp(as.matrix(df_alpha_covariates[, 2:n_covariates]) %*% cox_b1$coeff)

      cum_haz_whole_fulltime = cox_whole$cum_haz$cum_haz
      survivalwhole_cov = exp(as.matrix(df_alpha_covariates[, 2:n_covariates]) %*% cox_whole$coeff)
    }

    if(num_of_cores > 1){
      library(foreach)

      ## num_of_cores set-up
      cores = num_of_cores
      cl = snow::makeCluster(cores[1])
      doSNOW::registerDoSNOW(cl)
      pb = txtProgressBar(max = observed_t2, style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      opts = list(progress = progress)

      alpha_cov_list = foreach(counter_i = 1:observed_t2, .options.snow = opts, .combine = 'c', .packages = 'pracma')%dopar%{
        index_now_full = (counter_i - 1) * n_covariates + 1:n_covariates
        index_now_uncut = (counter_i - 1) * n_covariates + 1:n_covariates_cov
        index_now = (counter_i - 1) * n_covariates_cov + 1:n_covariates_cov

        alpha_cov_tmp = matrix(0, counter_i * n_covariates_cov, n_covariates_cov)
        alpha_cov_tmp[index_now, ] = alpha_var_tmp[index_now_uncut, ]

        ## fill off-diagonal
        if(counter_i > 1){
          sub_i = df_alpha_covariates[rank_t2[counter_i]:m, ]

          if(variance_method == "new"){
            survival1_cov_now = survival1_cov[rank_t2[counter_i]:m]
            survivalwhole_cov_now = survivalwhole_cov[rank_t2[counter_i]:m]
          }

          p2 = as.vector(1 / (1 + exp(-sub_i %*% alpha_mat[, counter_i])))
          sub_i_cov_2 = sub_i %*% alpha_var_tmp[index_now_full, ] # counter_i = cov_2

          index_now_j_full = 1:n_covariates
          index_now_j_uncut = 1:n_covariates_cov
          index_now_j = 1:n_covariates_cov
          if(rank_t2[counter_i] == m){
            for(counter_j in 1:(counter_i - 1)){
              if(variance_method == "new"){
                surv_n1_1 = exp(survival1_cov_now * (cum_haz_1_fulltime[counter_j] - cum_haz_1_fulltime[counter_i]))
                surv_n1_whole = exp(survivalwhole_cov_now * (cum_haz_whole_fulltime[counter_j] - cum_haz_whole_fulltime[counter_i]))
              }
              p1 = 1/(1 + exp(-as.vector(sub_i %*% alpha_mat[, counter_j])))

              if(variance_method == 'new'){
                prob = surv_n1_1/surv_n1_whole * p1 * (1 - p2)
              }else{
                prob = p1 * (1 - p2)
              }

              alpha_cov_tmp[index_now_j, ] = alpha_var[index_now_j_uncut, ] %*% (sub_i * prob) %*% sub_i_cov_2
              index_now_j = index_now_j + n_covariates_cov
              index_now_j_uncut = index_now_j_uncut + n_covariates
              index_now_j_full = index_now_j_full + n_covariates
            }
          }else{
            for(counter_j in 1:(counter_i - 1)){
              if(variance_method == "new"){
                surv_n1_1 = exp(survival1_cov_now * (cum_haz_1_fulltime[counter_j] - cum_haz_1_fulltime[counter_i]))
                surv_n1_whole = exp(survivalwhole_cov_now * (cum_haz_whole_fulltime[counter_j] - cum_haz_whole_fulltime[counter_i]))
              }
              p1 = 1/(1 + exp(-as.vector(sub_i %*% alpha_mat[, counter_j])))

              if(variance_method == 'new'){
                prob = surv_n1_1/surv_n1_whole * p1 * (1 - p2)
              }else{
                prob = p1 * (1 - p2)
              }

              alpha_cov_tmp[index_now_j, ] = tcrossprod(alpha_var[index_now_j_uncut, ], sub_i * prob) %*% sub_i_cov_2
              index_now_j = index_now_j + n_covariates_cov
              index_now_j_uncut = index_now_j_uncut + n_covariates
              index_now_j_full = index_now_j_full + n_covariates
            }
          }
        }
        gc()
        return(list(alpha_cov_tmp))
      }
      snow::stopCluster(cl)
      if(timer){pracma::fprintf('\n')}
      index_now = 1:n_covariates_cov
      for(counter_i in 1:observed_t2){
        alpha_cov[1:(counter_i * n_covariates_cov), index_now] = alpha_cov_list[[counter_i]]
        index_now = index_now + n_covariates_cov
      }
    }else{
      if(timer){
        pracma::fprintf('| variance estimation        30        40        50        60        70        80        90    100 |\n')
        loop_count = (1:observed_t2) * (1:observed_t2 - 1) / 2
        counter_total = observed_t2 * (observed_t2 - 1) / 2
        cum_bar_num = my_eva_fun(list(1:space, 1:space / space * counter_total), loop_count, rule = '0')
        bar_num = diff(c(0, cum_bar_num))
      }
      for(counter_i in 1:observed_t2){
        index_now = (counter_i - 1) * n_covariates_cov + 1:n_covariates_cov
        index_now_full = (counter_i - 1) * n_covariates + 1:n_covariates
        index_now_uncut = (counter_i - 1) * n_covariates + 1:n_covariates_cov

        ## fill diagonal block
        alpha_cov[index_now, index_now] = alpha_var_tmp[index_now_uncut, ]

        ## fill off-diagonal
        if(counter_i > 1){
          sub_i = df_alpha_covariates[rank_t2[counter_i]:m, ]
          if(variance_method == "new"){
            survival1_cov_now = survival1_cov[rank_t2[counter_i]:m]
            survivalwhole_cov_now = survivalwhole_cov[rank_t2[counter_i]:m]
          }

          p2 = as.vector(1 / (1 + exp(-sub_i %*% alpha_mat[, counter_i])))
          sub_i_cov_2 = sub_i %*% alpha_var_tmp[index_now_full, ] # counter_i = cov_2

          index_now_j_full = 1:n_covariates
          index_now_j_uncut = 1:n_covariates_cov
          index_now_j = 1:n_covariates_cov
          alpha_cov_tmp = matrix(0, (counter_i - 1) * n_covariates_cov, n_covariates_cov)
          if(rank_t2[counter_i] == m){
            for(counter_j in 1:(counter_i - 1)){
              if(variance_method == "new"){
                surv_n1_1 = exp(survival1_cov_now * (cum_haz_1_fulltime[counter_j] - cum_haz_1_fulltime[counter_i]))
                surv_n1_whole = exp(survivalwhole_cov_now * (cum_haz_whole_fulltime[counter_j] - cum_haz_whole_fulltime[counter_i]))
              }
              p1 = 1/(1 + exp(-as.vector(sub_i %*% alpha_mat[, counter_j])))

              if(variance_method == "new"){
                prob = surv_n1_1/surv_n1_whole * p1 * (1 - p2)
              }else{
                prob = p1 * (1 - p2)
              }

              alpha_cov_tmp[index_now_j, ] = alpha_var[index_now_j_uncut, ] %*% (sub_i * prob) %*% sub_i_cov_2
              index_now_j = index_now_j + n_covariates_cov
              index_now_j_uncut = index_now_j_uncut + n_covariates
              index_now_j_full = index_now_j_full + n_covariates
            }
          }else{
            for(counter_j in 1:(counter_i - 1)){
              if(variance_method == "new"){
                surv_n1_1 = exp(survival1_cov_now * (cum_haz_1_fulltime[counter_j] - cum_haz_1_fulltime[counter_i]))
                surv_n1_whole = exp(survivalwhole_cov_now * (cum_haz_whole_fulltime[counter_j] - cum_haz_whole_fulltime[counter_i]))
              }
              p1 = 1/(1 + exp(-as.vector(sub_i %*% alpha_mat[, counter_j])))

              if(variance_method == "new"){
                prob = surv_n1_1/surv_n1_whole * p1 * (1 - p2)
              }else{
                prob = p1 * (1 - p2)
              }

              alpha_cov_tmp[index_now_j, ] = tcrossprod(alpha_var[index_now_j_uncut, ], sub_i * prob) %*% sub_i_cov_2
              index_now_j = index_now_j + n_covariates_cov
              index_now_j_uncut = index_now_j_uncut + n_covariates
              index_now_j_full = index_now_j_full + n_covariates
            }
          }
          alpha_cov[1:(index_now[1] - 1), index_now] = alpha_cov_tmp
        }
        if(timer && bar_num[counter_i] > 0){for(i in 1:bar_num[counter_i]){pracma::fprintf('-')}}
      }
      if(timer){pracma::fprintf('\n')}
    }
    gdata::lowerTriangle(alpha_cov) = gdata::upperTriangle(alpha_cov, byrow = TRUE)
    return(list(time = unique_T2, coeff = alpha_mat, cov = alpha_cov, sick_alive = sick_alive, healthy_alive = healthy_alive, alive = alive, converged_alpha = converged_alpha, exact_time = exact_time))
  }else{
    return(list(time = unique_T2, coeff = alpha_mat, sick_alive = sick_alive, healthy_alive = healthy_alive, alive = alive, converged_alpha = converged_alpha, exact_time = exact_time))
  }
}
#' @export
mycoxph = function(time, observed, covariates, get_variance = TRUE){
  # time = time_b0; observed = observed_b0; covariates = covariates
  # time = time_b1; observed = observed_b1; covariates = as.matrix(covariates[df$d1, ])

  if(sum(observed) == 0){
    return(NULL)
  }else{
    Surv_object = survival::Surv(time[, 1], time[, 2], observed)
    if(var(covariates[, 1])==0){
      stop('All exposures are identical.')
    }
    cox = survival::coxph(Surv_object ~ covariates, method = "breslow", timefix = FALSE)
    coeff = cox$coefficients

    cum_haz = my_basehaz(time, observed, covariates, cox)

    # cum_haz = basehaz(cox, centered = FALSE)
    # colnames(cum_haz)[1] = 'cum_haz'
    # haz_logi = (diff(c(0, cum_haz$cum_haz)) != 0)
    # cum_haz = cum_haz[haz_logi, ]

    result = list()
    result$coeff = coeff
    result$cum_haz = cum_haz
    result$sub_df = data.frame(tstart = time[, 1], tend = time[, 2], observed = observed, covariates)
    AsymVariance = sum(c('a', 'A', 'asym', 'asymptotic', 'asymptotical', 'Asym', 'Asymptotic', 'Asymptotical') %in% get_variance) > 0

    if(AsymVariance){result$cov = inv_coxinformation(result$sub_df, coeff, cum_haz)}
  }
  return(result)
}

## asymptotics
#' @export
inv_coxinformation = function(df_, coeff_, cum_haz_){
  ## example: a
  # df_ = df_all$a
  # coeff_ = coeff$a
  # cum_haz_ = cum_haz$a

  ## example: b0
  # df_ = df_all$b0
  # coeff_ = coeff$b0
  # cum_haz_ = cum_haz$b0

  ## example: b1
  # df_ = df_all$b1
  # coeff_ = coeff$b1
  # cum_haz_ = cum_haz$b1

  coxinf = list()
  ## preprocess
  m = dim(df_)[1]
  n_col = dim(df_)[2]
  coeff_len = length(coeff_)
  covariates = df_[, 3 + 1:coeff_len]
  time_data = df_[, 1:3]

  ## upper-left
  if(coeff_len==1){
    exp_tmp = exp(covariates*coeff_)
  }else{
    covariates = as.matrix(covariates)
    exp_tmp = exp(covariates%*%coeff_)
  }
  tmp0 = exp_tmp*(my_eva_fun(cum_haz_, time_data[, 2]) - my_eva_fun(cum_haz_, time_data[, 1]))
  if(coeff_len == 1){
    coxinf$upper_left = as.matrix(sum(covariates^2*tmp0))
  }else{
    coxinf$upper_left = t(covariates*pracma::repmat(tmp0, 1, coeff_len))%*%as.matrix(covariates)
  }

  ## lower-left ## 1 for start; 2 for end;
  rank_time_2 = rank(time_data[time_data[, 3] == 1, 2], ties.method = 'min')
  sort_time_2 = sort(time_data[, 2], index.return = TRUE)
  sorted_time_df_2 = time_data[sort_time_2$ix, ]

  if(coeff_len==1){
    sorted_covariates = covariates[sort_time_2$ix]
    cov_exp_tmp = sorted_covariates*as.vector(exp(sorted_covariates*coeff_))
    tmp_lower_left = (sum(cov_exp_tmp) - c(0, cumsum(cov_exp_tmp)))[-(m+1)]
    tmp_lower_left = tmp_lower_left[sorted_time_df_2[, 3] == 1]
    tmp2 = tmp_lower_left[rank_time_2]
  }else{
    sorted_covariates = covariates[sort_time_2$ix,]
    cov_exp_tmp = sorted_covariates*pracma::repmat(exp(sorted_covariates%*%coeff_), 1, coeff_len)
    tmp_lower_left = apply(cov_exp_tmp, 2, function(x) sum(x) - c(0, cumsum(x)))[-(m+1), ]
    tmp_lower_left = tmp_lower_left[sorted_time_df_2[, 3] == 1, ]
    tmp2 = tmp_lower_left[rank_time_2, ]
  }

  sort_time_1 = sort(time_data[, 1], index.return = TRUE)
  sorted_time_df_1 = time_data[sort_time_1$ix, ]
  important_index = get_position(df_$tend[df_$observed==1], sort_time_1$x)

  if(coeff_len==1){
    sorted_covariates = covariates[sort_time_1$ix]
    cov_exp_tmp = sorted_covariates*as.vector(exp(sorted_covariates*coeff_))
    tmp_lower_left = (sum(cov_exp_tmp) - c(0, cumsum(cov_exp_tmp))) # use important index
    tmp1 = tmp_lower_left[important_index]
  }else{
    sorted_covariates = covariates[sort_time_1$ix,]
    cov_exp_tmp = sorted_covariates*pracma::repmat(exp(sorted_covariates%*%coeff_), 1, coeff_len)
    tmp_lower_left = apply(cov_exp_tmp, 2, function(x) sum(x) - c(0, cumsum(x)))
    tmp1 = tmp_lower_left[important_index, ]
  }
  coxinf$lower_left = tmp2 - tmp1

  ## inverse of lower-right
  tmp_inv_lower_right = diff(c(0, cum_haz_$cum_haz))^2
  tmp_inv_lower_right = tmp_inv_lower_right[tmp_inv_lower_right!=0]
  coxinf$inv_lower_right = tmp_inv_lower_right[rank_time_2]

  ## inverse of coxinf
  # quadratic_core = (A-BD^{-1}C)^{-1}
  # linear_core = BD^{-1} vector or matrix
  linear_core = as.matrix(t(coxinf$lower_left * coxinf$inv_lower_right))
  # quadratic_core = cov_
  quadratic_core = as.matrix(solve(coxinf$upper_left - linear_core %*% coxinf$lower_left))

  inv_coxinf = list()
  inv_coxinf$upper_left = quadratic_core
  inv_coxinf$upper_right = -quadratic_core %*% linear_core
  inv_coxinf$lower_left = t(inv_coxinf$upper_right)
  inv_coxinf$inv_coxinf_lower_right = coxinf$inv_lower_right

  inv_coxinf$linear_core = t(linear_core)
  inv_coxinf$t_linear_core = linear_core
  inv_coxinf$quadratic_core = quadratic_core
  inv_coxinf$coxinf = coxinf
  return(inv_coxinf)
}
#' @export
compute_variance = function(get_DE, get_IE, intervention, cal_level, estimation_alpha, cox_b0, cox_b1, b0_time, b1_time){
  if(get_DE){za_iv1 = intervention[1]; zb_iv1 = intervention[2]; za_iv2 = intervention[2]; zb_iv2 = intervention[2];}
  if(get_IE){za_iv1 = intervention[1]; zb_iv1 = intervention[1]; za_iv2 = intervention[1]; zb_iv2 = intervention[2];}

  # 1 and 2 in pd1 and pd2 stand for two different types of intervention.
  pd_iv1 = get_pd(za_iv1, zb_iv1, estimation_alpha, cal_level, cox_b0, cox_b1, b0_time, b1_time) # za = za0
  pd_iv2 = get_pd(za_iv2, zb_iv2, estimation_alpha, cal_level, cox_b0, cox_b1, b0_time, b1_time) # zb = zb0

  ## alpha variance
  pd_alpha = pd_iv1$alpha_tmp - pd_iv2$alpha_tmp
  alpha_variance = get_alpha_variance(pd_alpha, estimation_alpha$cov)

  ## beta variance
  if(!is.null(cox_b0)){
    pd_beta_vec_0 = pd_iv1$beta_0 - pd_iv2$beta_0
    pd_cum_beta_0 = pd_iv1$beta_baseline_0 - pd_iv2$beta_baseline_0
    beta_variance_0 = list(get_beta_variance(pd_beta_vec_0, pd_cum_beta_0, cox_b0$cov), cox_b0$cum_haz$time)
  }else{beta_variance_0 = list(rep(0, length(alpha_variance)), estimation_alpha$time)}

  if(!is.null(cox_b1)){
    pd_beta_vec_1 = pd_iv1$beta_1 - pd_iv2$beta_1
    pd_cum_beta_1 = pd_iv1$beta_baseline_1 - pd_iv2$beta_baseline_1
    beta_variance_1 = list(get_beta_variance(pd_beta_vec_1, pd_cum_beta_1, cox_b1$cov), cox_b1$cum_haz$time)
  }else{beta_variance_1 = list(rep(0, length(alpha_variance)), estimation_alpha$time)}

  beta_variance = my_eva_fun(beta_variance_0, estimation_alpha$time, rule = '0') + my_eva_fun(beta_variance_1, estimation_alpha$time, rule = '0')

  variance = alpha_variance + beta_variance
  return(data.frame(variance = variance, alpha_variance = alpha_variance, beta_variance = beta_variance))
}
#' @export
get_pd = function(za_iv, zb_iv, estimation_alpha, cal_level, cox_b0, cox_b1, b0_time, b1_time){
  # za_iv = za_iv1; zb_iv = zb_iv1
  # za_iv = za_iv2; zb_iv = zb_iv2

  pd = list()
  zaX = c(za_iv, cal_level); exp_zaX0 = exp(sum(zaX * cox_b0$coeff)); exp_zaX1 = exp(sum(zaX * cox_b1$coeff))
  zbX = c(1, zb_iv, cal_level[cal_level != 0]); tmp_zbX = as.vector(zbX %*% estimation_alpha$coeff[1:length(zbX), ])
  if(length(estimation_alpha$exact_time$time) == length(estimation_alpha$time)){
    full_zbX = tmp_zbX
  }else{
    full_zbX = approx(x = estimation_alpha$time, y = tmp_zbX, xout = estimation_alpha$exact_time$time, method = 'linear', rule = 2)$y
  }
  exp_zbX  = exp(full_zbX)
  beta_denominator = 1/(1 + exp_zbX)

  group = as.numeric(estimation_alpha$exact_time$time %in% cox_b1$cum_haz$time)
  small_lambda_n1_0 = diff(c(0, cox_b0$cum_haz$cum_haz))
  small_lambda_n1_1 = diff(c(0, cox_b1$cum_haz$cum_haz))

  ## beta
  # baseline -- n1 = 0
  pd$beta_baseline_0 = beta_denominator[group == 0] * exp_zaX0
  # baseline -- n1 = 1
  pd$beta_baseline_1 = (1 - beta_denominator[group == 1]) * exp_zaX1

  # covariates -- n1 = 0
  pd$beta_0 = tcrossprod(cumsum(pd$beta_baseline_0 * small_lambda_n1_0), zaX)
  # covariates -- n1 = 1
  pd$beta_1 = tcrossprod(cumsum(pd$beta_baseline_1 * small_lambda_n1_1), zaX)

  ## alpha
  # alpha_core = beta_denominator_alpha * (1 - beta_denominator_alpha)
  alpha_core = beta_denominator * (1 - beta_denominator)
  alpha_tmp = matrix(0, length(group), length(zbX))
  # n1 = 0
  alpha_tmp[group == 0, ] = -tcrossprod(alpha_core[group == 0] * small_lambda_n1_0 * exp_zaX0, zbX)
  # n1 = 1
  alpha_tmp[group == 1, ] = tcrossprod(alpha_core[group == 1] * small_lambda_n1_1 * exp_zaX1, zbX)

  # "true" alpha
  pd$alpha_tmp = matrix(0, length(estimation_alpha$exact_time$rank_exact_t2), length(zbX))
  for(i in 1:(dim(pd$alpha_tmp)[1] - 1)){
    index_now = estimation_alpha$exact_time$rank_exact_t2[c(i, i+1)]
    if(index_now[2] - index_now[1] == 1){
      pd$alpha_tmp[i, ] = alpha_tmp[index_now[1], ]
    }else{
      pd$alpha_tmp[i, ] = colSums(alpha_tmp[index_now[1] : (index_now[2] - 1), ])
    }
  }

  return(pd)
}
#' @export
get_alpha_variance = function(pd_alpha, alpha_cov){
  # alpha_cov = estimation_alpha$cov
  n_row = dim(pd_alpha)[1]
  pd_alpha_not0 = (colSums(pd_alpha) != 0)
  cut_it = sum(pd_alpha_not0) == 0
  n_covariates_uncut = dim(alpha_cov)[1] / dim(pd_alpha)[1]
  if(cut_it){
    pd_alpha = as.matrix(pd_alpha[, pd_alpha_not0])
    chosed_index_tmp = which(pd_alpha_not0)
  }
  n_covariates = dim(pd_alpha)[2]

  pd_alpha_mat = matrix(0, n_row, n_row * n_covariates)
  pd_alpha_tmp = rep(0, n_row * n_covariates)
  chosed_index = rep(0, n_row * n_covariates)
  now_index = 1:n_covariates

  for(i in 1:n_row){
    pd_alpha_tmp[now_index] = pd_alpha[i, ]
    pd_alpha_mat[i, ] = pd_alpha_tmp

    if(cut_it){
      chosed_index[now_index] = chosed_index_tmp
      chosed_index_tmp = chosed_index_tmp + n_covariates_uncut
    }
    now_index = now_index + n_covariates
  }
  if(cut_it){alpha_cov = alpha_cov[chosed_index, chosed_index]}
  alpha_variance = rowSums(pd_alpha_mat %*% alpha_cov * pd_alpha_mat)
  alpha_variance[alpha_variance < 0] = 0
  return(alpha_variance)
}
#' @export
get_beta_variance = function(pd_beta_vec, pd_cum_beta, beta_cov){
  # pd_beta_vec = pd_beta_vec_0; pd_cum_beta = pd_cum_beta_0; beta_cov = cox_b0$cov;
  # pd_beta_vec = pd_beta_vec_1; pd_cum_beta = pd_cum_beta_1; beta_cov = cox_b1$cov;

  n_col_b = sqrt(length(beta_cov$upper_left))
  matrix_pd_cum_beta = cbind(pd_cum_beta)[, rep(1, n_col_b)]
  beta_tmp_1 = pd_beta_vec%*%beta_cov$upper_left
  beta_tmp_2 = apply(beta_cov$lower_left*matrix_pd_cum_beta, 2, cumsum)
  beta_part1 = rowSums(pd_beta_vec*(beta_tmp_1 + 2*beta_tmp_2))

  beta_tmp_3 = cumsum(pd_cum_beta^2*beta_cov$inv_coxinf_lower_right)
  beta_tmp_4 = rowSums(beta_tmp_2 * apply(beta_cov$linear_core*matrix_pd_cum_beta, 2, cumsum))
  beta_part2 = beta_tmp_3 - beta_tmp_4
  beta_variance = beta_part1 + beta_part2
  return(beta_variance)
}

## sensitivity analysis
#' @export
do_sen_ana = function(get_DE, get_IE, intervention, cal_level, estimation_alpha, small_cox_b0, small_cox_b1){
  if(get_DE){za_iv1 = intervention[1]; zb_iv1 = intervention[2]; za_iv2 = intervention[2]; zb_iv2 = intervention[2];}
  if(get_IE){za_iv1 = intervention[1]; zb_iv1 = intervention[1]; za_iv2 = intervention[1]; zb_iv2 = intervention[2];}
  alpha_hat = list(time = estimation_alpha$time, coeff = estimation_alpha$coeff)
  cox_hat_0 = list(coeff = small_cox_b0$coeff, cum_haz = small_cox_b0$cum_haz)
  cox_hat_1 = list(coeff = small_cox_b1$coeff, cum_haz = small_cox_b1$cum_haz)

  ## case 1
  counter = 0
  gamma_list = seq(-1, 1, by = 0.2)
  case1 = list(gamma = gamma_list, effect = vector(mode = 'list', length = length(gamma_list)))
  for(gamma in gamma_list){
    counter = counter + 1

    alpha_U = gamma
    beta_U = gamma
    gamma_Z = abs(gamma)
    gamma_Z2 = abs(gamma)
    gamma_n1 = abs(gamma)

    ## adjust alpha
    alpha_hat$coeff[1, ] = estimation_alpha$coeff[1, ] - 1/2 * alpha_U^2
    alpha_hat$coeff[2, ] = estimation_alpha$coeff[2, ] - alpha_U * gamma_Z2

    ## adjust beta
    cox_hat_0$coeff = small_cox_b0$coeff - beta_U * gamma_Z
    cox_hat_1$coeff = small_cox_b1$coeff - beta_U * gamma_Z

    cox_hat_0$cum_haz$cum_haz = small_cox_b0$cum_haz$cum_haz/exp(beta_U ^ 2/2)
    cox_hat_1$cum_haz$cum_haz = small_cox_b1$cum_haz$cum_haz/exp(beta_U ^ 2/2)

    ## adjust counterfactual hazard
    counterfactual_hazard_iv1 = get_counterfactual_hazard(za_iv1, zb_iv1, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard_iv2 = get_counterfactual_hazard(za_iv2, zb_iv2, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard = counterfactual_hazard_iv1 - counterfactual_hazard_iv2

    case1$effect[[counter]] = data.frame(time = estimation_alpha$time, effect = counterfactual_hazard)
  }

  ## case 2
  counter = 0
  gamma_list = seq(0, 1, by = 0.1)
  case2 = list(gamma = gamma_list, effect = vector(mode = 'list', length = length(gamma_list)))
  for(gamma in gamma_list){
    counter = counter + 1

    alpha_U = log(1.5)
    beta_U = log(1.5)
    gamma_Z = abs(gamma)
    gamma_Z2 = abs(gamma)
    gamma_n1 = abs(gamma)

    ## adjust alpha
    alpha_hat$coeff[1, ] = estimation_alpha$coeff[1, ] - 1/2 * alpha_U^2
    alpha_hat$coeff[2, ] = estimation_alpha$coeff[2, ] - alpha_U * gamma_Z2

    ## adjust beta
    cox_hat_0$coeff = small_cox_b0$coeff - beta_U * gamma_Z
    cox_hat_1$coeff = small_cox_b1$coeff - beta_U * gamma_Z

    cox_hat_0$cum_haz$cum_haz = small_cox_b0$cum_haz$cum_haz/exp(beta_U ^ 2/2)
    cox_hat_1$cum_haz$cum_haz = small_cox_b1$cum_haz$cum_haz/exp(beta_U ^ 2/2)

    ## adjust counterfactual hazard
    counterfactual_hazard_iv1 = get_counterfactual_hazard(za_iv1, zb_iv1, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard_iv2 = get_counterfactual_hazard(za_iv2, zb_iv2, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard = counterfactual_hazard_iv1 - counterfactual_hazard_iv2

    case2$effect[[counter]] = data.frame(time = estimation_alpha$time, effect = counterfactual_hazard)
  }

  ## case 3
  counter = 0
  gamma_list = seq(0, 1, by = 0.1)
  case3 = list(gamma = gamma_list, effect = vector(mode = 'list', length = length(gamma_list)))
  for(gamma in gamma_list){
    counter = counter + 1

    alpha_U = -log(1.5)
    beta_U = -log(1.5)
    gamma_Z = abs(gamma)
    gamma_Z2 = abs(gamma)
    gamma_n1 = abs(gamma)

    ## adjust alpha
    alpha_hat$coeff[1, ] = estimation_alpha$coeff[1, ] - 1/2 * alpha_U^2
    alpha_hat$coeff[2, ] = estimation_alpha$coeff[2, ] - alpha_U * gamma_Z2

    ## adjust beta
    cox_hat_0$coeff = small_cox_b0$coeff - beta_U * gamma_Z
    cox_hat_1$coeff = small_cox_b1$coeff - beta_U * gamma_Z

    cox_hat_0$cum_haz$cum_haz = small_cox_b0$cum_haz$cum_haz/exp(beta_U ^ 2/2)
    cox_hat_1$cum_haz$cum_haz = small_cox_b1$cum_haz$cum_haz/exp(beta_U ^ 2/2)

    ## adjust counterfactual hazard
    counterfactual_hazard_iv1 = get_counterfactual_hazard(za_iv1, zb_iv1, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard_iv2 = get_counterfactual_hazard(za_iv2, zb_iv2, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard = counterfactual_hazard_iv1 - counterfactual_hazard_iv2

    case3$effect[[counter]] = data.frame(time = estimation_alpha$time, effect = counterfactual_hazard)
  }

  ## case 4
  counter = 0
  gamma_list = seq(-1, 1, by = 0.2)
  case4 = list(gamma = gamma_list, effect = vector(mode = 'list', length = length(gamma_list)))
  for(gamma in gamma_list){
    counter = counter + 1

    alpha_U = gamma
    beta_U = gamma
    gamma_Z = 1.2
    gamma_Z2 = 1.2
    gamma_n1 = 1.2

    ## adjust alpha
    alpha_hat$coeff[1, ] = estimation_alpha$coeff[1, ] - 1/2 * alpha_U^2
    alpha_hat$coeff[2, ] = estimation_alpha$coeff[2, ] - alpha_U * gamma_Z2

    ## adjust beta
    cox_hat_0$coeff = small_cox_b0$coeff - beta_U * gamma_Z
    cox_hat_1$coeff = small_cox_b1$coeff - beta_U * gamma_Z

    cox_hat_0$cum_haz$cum_haz = small_cox_b0$cum_haz$cum_haz/exp(beta_U ^ 2/2)
    cox_hat_1$cum_haz$cum_haz = small_cox_b1$cum_haz$cum_haz/exp(beta_U ^ 2/2)

    ## adjust counterfactual hazard
    counterfactual_hazard_iv1 = get_counterfactual_hazard(za_iv1, zb_iv1, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard_iv2 = get_counterfactual_hazard(za_iv2, zb_iv2, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard = counterfactual_hazard_iv1 - counterfactual_hazard_iv2

    case4$effect[[counter]] = data.frame(time = estimation_alpha$time, effect = counterfactual_hazard)
  }

  ## case 5
  counter = 0
  gamma_list = seq(-1, 1, by = 0.2)
  case5 = list(gamma = gamma_list, effect = vector(mode = 'list', length = length(gamma_list)))
  for(gamma in gamma_list){
    counter = counter + 1

    alpha_U = 0
    beta_U = 0
    gamma_Z = gamma
    gamma_Z2 = gamma
    gamma_n1 = gamma

    ## adjust alpha
    alpha_hat$coeff[1, ] = estimation_alpha$coeff[1, ] - 1/2 * alpha_U^2
    alpha_hat$coeff[2, ] = estimation_alpha$coeff[2, ] - alpha_U * gamma_Z2

    ## adjust beta
    cox_hat_0$coeff = small_cox_b0$coeff - beta_U * gamma_Z
    cox_hat_1$coeff = small_cox_b1$coeff - beta_U * gamma_Z

    cox_hat_0$cum_haz$cum_haz = small_cox_b0$cum_haz$cum_haz/exp(beta_U ^ 2/2)
    cox_hat_1$cum_haz$cum_haz = small_cox_b1$cum_haz$cum_haz/exp(beta_U ^ 2/2)

    ## adjust counterfactual hazard
    counterfactual_hazard_iv1 = get_counterfactual_hazard(za_iv1, zb_iv1, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard_iv2 = get_counterfactual_hazard(za_iv2, zb_iv2, cal_level, alpha_hat, cox_hat_0, cox_hat_1)
    counterfactual_hazard = counterfactual_hazard_iv1 - counterfactual_hazard_iv2

    case5$effect[[counter]] = data.frame(time = estimation_alpha$time, effect = counterfactual_hazard)
  }
  return(list(case1 = case1, case2 = case2, case3 = case3, case4 = case4, case5 = case5))
}

## counterfactual hazard
#' @export
get_counterfactual_hazard = function(za_iv, zb_iv, cal_level, estimation_alpha, cox_b0, cox_b1){
  # za_iv = za_iv1; zb_iv = zb_iv1;
  # za_iv = za_iv2; zb_iv = zb_iv2;
  # cox_b0 = small_cox_b0; cox_b1 = small_cox_b1

  intercept = 1
  w_prob = as.vector(1/(1 + exp(-crossprod(c(intercept, zb_iv, cal_level), estimation_alpha$coeff))))

  group_0_time = estimation_alpha$time %in% cox_b0$cum_haz$time
  group_1_time = estimation_alpha$time %in% cox_b1$cum_haz$time

  if(sum(group_0_time) > 0){
    w0 = 1 - w_prob[group_0_time]
    dLbase_0 = diff(c(0, cox_b0$cum_haz$cum_haz))
    dL0 = dLbase_0 * exp(sum(c(za_iv, cal_level) * cox_b0$coeff))
    n1_0 = my_eva_fun(list(cumsum(w0 * dL0), cox_b0$cum_haz$time), estimation_alpha$time, rule = '0')
  }else{
    n1_0 = rep(0, length(group_0_time))
  }

  if(sum(group_1_time) > 0){
    w1 = w_prob[group_1_time]
    dLbase_1 = diff(c(0, cox_b1$cum_haz$cum_haz))
    dL1 = dLbase_1 * exp(sum(c(za_iv, cal_level) * cox_b1$coeff))
    n1_1 = my_eva_fun(list(cumsum(w1 * dL1), cox_b1$cum_haz$time), estimation_alpha$time, rule = '0')
  }else{
    n1_1 = rep(0, length(group_1_time))
  }

  counterfactual_hazard = n1_0 + n1_1
  return(counterfactual_hazard)
}
#' @export
estimate_effect = function(df, effect, intervention, cal_level, sen_ana, get_variance, boot_times, timer, num_of_cores, unique_T2, b0_time, b1_time, variance_method, threshold, HO = FALSE){
  ## beta part
  # auxiliary
  m = dim(df)[1]
  n_col = dim(df)[2]
  num_covariates = n_col - 4
  covariates = as.matrix(df[, 4+1:num_covariates])

  # b0
  time_b0 = cbind(rep(0, m), df$T1)
  observed_b0 = df$d2 & (df$d1 == FALSE)
  cox_b0 = mycoxph(time_b0, observed_b0, covariates, get_variance)

  # b1
  df_b1 = df[df$d1, ]
  time_b1 = df_b1[, c(1, 2)]
  observed_b1 = df_b1$d2
  cox_b1 = mycoxph(time_b1, observed_b1, as.matrix(covariates[df$d1, ]), get_variance)

  # whole data
  time_whole = df$T2
  cox_whole = mycoxph(cbind(rep(0, m), df$T2), df$d2, as.matrix(covariates), get_variance = NULL)

  # make it small
  #-----------------------------------------------------------------------------------------------#
  # small cox will be used while computing the covariance of alphas and the counterfactual hazard #
  # cox will be used while computing the variance                                                 #
  #-----------------------------------------------------------------------------------------------#

  small_cox_b0 = make_small(cox_b0, b0_time)
  small_cox_b1 = make_small(cox_b1, b1_time)
  small_cox_whole = make_small(cox_whole, unique_T2)

  ## alpha part
  estimation_alpha = estimate_alpha(df, cal_level, small_cox_b1, small_cox_whole, unique_T2, get_variance, timer, num_of_cores, variance_method, threshold)

  ## get counterfactual hazard
  AsymVariance = sum(c('a', 'A', 'asym', 'asymptotic', 'asymptotical', 'Asym', 'Asymptotic', 'Asymptotical') %in% get_variance) > 0

  get_DE = sum(c('d', 'D', 'de', 'De', 'DE', 'direct effect', 'Direct effect', 'Direct Effect') %in% effect) > 0
  get_IE = sum(c('i', 'I', 'ie', 'Ie', 'IE', 'indirect effect', 'Indirect effect', 'Indirect Effect') %in% effect) > 0
  result = list()

  ## Direct effect
  if(get_DE){
    za_iv1 = intervention[1]; zb_iv1 = intervention[2]; za_iv2 = intervention[2]; zb_iv2 = intervention[2];
    counterfactual_hazard_iv1 = get_counterfactual_hazard(za_iv1, zb_iv1, cal_level, estimation_alpha, small_cox_b0, small_cox_b1)
    counterfactual_hazard_iv2 = get_counterfactual_hazard(za_iv2, zb_iv2, cal_level, estimation_alpha, small_cox_b0, small_cox_b1)
    counterfactual_hazard = counterfactual_hazard_iv1 - counterfactual_hazard_iv2

    result$DE$effect = counterfactual_hazard
    result$DE$time = estimation_alpha$time
    if(AsymVariance){
      result$DE$variance = compute_variance(get_DE = TRUE, get_IE = FALSE, intervention, cal_level, estimation_alpha, cox_b0, cox_b1, b0_time, b1_time)
      result$DE$asym_lower = result$DE$effect - 1.96 * sqrt(result$DE$variance$variance)
      result$DE$asym_upper = result$DE$effect + 1.96 * sqrt(result$DE$variance$variance)
    }
    if(sen_ana){result$DE$sensitivity_analysis = do_sen_ana(get_DE = TRUE, get_IE = FALSE, intervention, cal_level, estimation_alpha, small_cox_b0, small_cox_b1)}
  }

  ## Indirect effect
  if(get_IE){
    za_iv1 = intervention[1]; zb_iv1 = intervention[1]; za_iv2 = intervention[1]; zb_iv2 = intervention[2];
    counterfactual_hazard_iv1 = get_counterfactual_hazard(za_iv1, zb_iv1, cal_level, estimation_alpha, small_cox_b0, small_cox_b1)
    counterfactual_hazard_iv2 = get_counterfactual_hazard(za_iv2, zb_iv2, cal_level, estimation_alpha, small_cox_b0, small_cox_b1)
    counterfactual_hazard = counterfactual_hazard_iv1 - counterfactual_hazard_iv2

    result$IE$effect = counterfactual_hazard
    result$IE$time = estimation_alpha$time

    if(HO){
      integrand = diff(c(0, counterfactual_hazard))
      weight = sqrt(estimation_alpha$sick_alive$number * estimation_alpha$healthy_alive$number)/estimation_alpha$alive$number
      result$IE$HO = data.frame(weight = weight, integrand = integrand, time = estimation_alpha$sick_alive$time)
    }

    if(AsymVariance){
      result$IE$variance = compute_variance(get_DE = FALSE, get_IE = TRUE, intervention, cal_level, estimation_alpha, cox_b0, cox_b1)
      result$IE$asym_lower = result$IE$effect - 1.96 * sqrt(result$IE$variance$variance)
      result$IE$asym_upper = result$IE$effect + 1.96 * sqrt(result$IE$variance$variance)
    }

    if(sen_ana){result$IE$sensitivity_analysis = do_sen_ana(get_DE = FALSE, get_IE = TRUE, intervention, cal_level, estimation_alpha, small_cox_b0, small_cox_b1)}
  }

  if(HO){
    result$variance$MA = estimation_alpha$cov
    if(length(cox_b0) > 1){result$variance$MB0 = solve(form_matrix(cox_b0$cov$upper_left, cox_b0$cov$lower_left, 1/cox_b0$cov$coxinf$inv_lower_right)/m)}
    if(length(cox_b1) > 1){result$variance$MB1 = solve(form_matrix(cox_b1$cov$upper_left, cox_b1$cov$lower_left, 1/cox_b1$cov$coxinf$inv_lower_right)/m)}
  }

  result$alive = estimation_alpha$alive
  result$healthy_alive = estimation_alpha$healthy_alive
  result$sick_alive = estimation_alpha$sick_alive
  result$converged_alpha = estimation_alpha$converged_alpha
  result$alpha = estimation_alpha$coeff
  result$cox_b0 = list(coeff = cox_b0$coeff, cum_haz = small_cox_b0$cum_haz)
  result$cox_b1 = list(coeff = cox_b1$coeff, cum_haz = small_cox_b1$cum_haz)
  return(result)
}

## plot function
plot_poly = function(y1, y2, x, color, density = NULL, angle = NULL){
  # y1 = df_asymp_DE$lower; y2 = df_asymp_DE$upper; x = df_asymp_DE$time; color = "dodgerblue"
  yy1 = c(rep(y1, each = 2)); yy1 = yy1[1:(length(yy1) - 1)]
  yy2 = c(rep(y2, each = 2)); yy2 = rev(yy2[1:(length(yy2)) - 1])
  xx = rep(x, each = 2); xx = xx[2:length(xx)]
  polygon(y = c(yy1, yy2), x = c(xx, rev(xx)), density = density, col = adjustcolor(color, alpha.f = 0.3), border = color, angle = angle)
}
plot_CHH2020 = function(result){
  result$IE$time = result$IE$time / 365.25
  result$DE$time = result$DE$time / 365.25
  if(!is.null(result$DE$boot_upper)){
    ylim_cumh_upper = max(result$IE$boot_upper, result$DE$boot_upper, result$IE$asym_upper, result$DE$asym_upper)
    ylim_cumh_lower = min(result$IE$boot_lower, result$DE$boot_lower, result$IE$asym_lower, result$DE$asym_lower)
  }else{
    ylim_cumh_upper = max(result$IE$asym_upper, result$DE$asym_upper)
    ylim_cumh_lower = min(result$IE$asym_lower, result$DE$asym_lower)
  }

  ylim_surv_lower = exp(-ylim_cumh_upper)
  ylim_surv_upper = exp(-ylim_cumh_lower)

  ## plot default
  cex.lab = 1.1
  cex.main = 1.1
  cex.axis = 1
  las = 1
  xlab = 'Time (years)'
  ylab_rho = expression(paste(rho[DE](t), ',', rho[IE](t)))
  ylab_Del = expression(paste(Delta[DE](t), ',', Delta[IE](t)))

  if(!is.null(result$DE$boot_upper)){
    ## bootstrap, surv
    df_asymp_IE = data.frame(cumhaz = exp(-result$IE$effect), time = result$IE$time, upper = exp(-result$IE$boot_upper), lower = exp(-result$IE$boot_lower))
    df_asymp_DE = data.frame(cumhaz = exp(-result$DE$effect), time = result$DE$time, upper = exp(-result$DE$boot_upper), lower = exp(-result$DE$boot_lower))

    plot(cumhaz ~ time, data = df_asymp_IE, type = "s", lwd = 2, ylim = c(ylim_surv_lower, ylim_surv_upper), col = adjustcolor("orange", alpha.f = 0.50), main = "Survival probability ratio \n Bootstrap CI", xlab = xlab, ylab = ylab_rho, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, las = las)
    lines(cumhaz ~ time, data = df_asymp_DE, type = "s", lwd = 2, col = adjustcolor("dodgerblue", alpha.f = 0.50))
    legend("bottomleft", legend = c("Indirect effect", "Direct effect"), fill = c(adjustcolor("orange", alpha.f = 0.10), adjustcolor("dodgerblue", alpha.f = 0.10)), bty = "n", border = c(adjustcolor("orange", alpha.f = 10), adjustcolor("dodgerblue", alpha.f = 10)))
    plot_poly(df_asymp_DE$lower, df_asymp_DE$upper, df_asymp_DE$time, "dodgerblue", NULL)
    plot_poly(df_asymp_IE$lower, df_asymp_IE$upper, df_asymp_IE$time, "orange", NULL)
    abline(h = 1, col = "grey")

    ## bootstrap, hazard
    df_asymp_IE = data.frame(cumhaz = result$IE$effect, time = result$IE$time, upper = result$IE$boot_upper, lower = result$IE$boot_lower)
    df_asymp_DE = data.frame(cumhaz = result$DE$effect, time = result$DE$time, upper = result$DE$boot_upper, lower = result$DE$boot_lower)

    plot(cumhaz ~ time, data = df_asymp_IE, type = "s", lwd = 2, ylim = c(ylim_cumh_lower, ylim_cumh_upper), col = adjustcolor("orange", alpha.f = 0.50), main = "Cumulative hazard difference \n Bootstrap CI",  xlab = xlab, ylab = ylab_Del, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, las = las)
    lines(cumhaz ~ time, data = df_asymp_DE, type = "s", lwd = 2, col = adjustcolor("dodgerblue", alpha.f = 0.50))
    legend("topleft", legend = c("Indirect effect", "Direct effect"), cex = 0.85, fill = c(adjustcolor("orange", alpha.f = 0.10), adjustcolor("dodgerblue", alpha.f = 0.10)), bty = "n", border = c(adjustcolor("orange", alpha.f = 10), adjustcolor("dodgerblue", alpha.f = 10)))
    plot_poly(df_asymp_DE$lower, df_asymp_DE$upper, df_asymp_DE$time, "dodgerblue", NULL)
    plot_poly(df_asymp_IE$lower, df_asymp_IE$upper, df_asymp_IE$time, "orange", NULL)
    abline(h = 0, col = "grey")
  }

  ## asymptotic, surv
  df_asymp_IE = data.frame(cumhaz = exp(-result$IE$effect), time = result$IE$time, upper = exp(-result$IE$asym_upper), lower = exp(-result$IE$asym_lower))
  df_asymp_DE = data.frame(cumhaz = exp(-result$DE$effect), time = result$DE$time, upper = exp(-result$DE$asym_upper), lower = exp(-result$DE$asym_lower))

  plot(cumhaz ~ time, data = df_asymp_IE, type = "s", lwd = 2, ylim = c(ylim_surv_lower, ylim_surv_upper), col = adjustcolor("orange", alpha.f = 0.50), main = "Survival probability ratio \n Asymptotic CI",  xlab = xlab, ylab = ylab_rho, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, las = las)
  lines(cumhaz ~ time, data = df_asymp_DE, type = "s", lwd = 2, col = adjustcolor("dodgerblue", alpha.f = 0.50))
  legend("bottomleft", legend = c("Indirect effect", "Direct effect"), cex = 0.85, fill = c(adjustcolor("orange", alpha.f = 0.10), adjustcolor("dodgerblue", alpha.f = 0.10)), bty = "n", border = c(adjustcolor("orange", alpha.f = 10), adjustcolor("dodgerblue", alpha.f = 10)))
  plot_poly(df_asymp_DE$lower, df_asymp_DE$upper, df_asymp_DE$time, "dodgerblue", NULL)
  plot_poly(df_asymp_IE$lower, df_asymp_IE$upper, df_asymp_IE$time, "orange", NULL)
  abline(h = 1, col = "grey")


  ## asymptotic, hazard
  df_asymp_IE = data.frame(cumhaz = result$IE$effect, time = result$IE$time, upper = result$IE$asym_upper, lower = result$IE$asym_lower)
  df_asymp_DE = data.frame(cumhaz = result$DE$effect, time = result$DE$time, upper = result$DE$asym_upper, lower = result$DE$asym_lower)

  plot(cumhaz ~ time, data = df_asymp_IE, type = "s", lwd = 2, ylim = c(ylim_cumh_lower, ylim_cumh_upper), col = adjustcolor("orange", alpha.f = 0.50), main = "Cumulative hazard difference \n Asymptotic CI", xlab = xlab, ylab = ylab_Del, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, las = las)
  lines(cumhaz ~ time, data = df_asymp_DE, type = "s", lwd = 2, col = adjustcolor("dodgerblue", alpha.f = 0.50))
  legend("topleft", legend = c("Indirect effect", "Direct effect"), cex = 0.85, fill = c(adjustcolor("orange", alpha.f = 0.10), adjustcolor("dodgerblue", alpha.f = 0.10)), bty = "n", border = c(adjustcolor("orange", alpha.f = 10), adjustcolor("dodgerblue", alpha.f = 10)))
  plot_poly(df_asymp_DE$lower, df_asymp_DE$upper, df_asymp_DE$time, "dodgerblue", NULL)
  plot_poly(df_asymp_IE$lower, df_asymp_IE$upper, df_asymp_IE$time, "orange", NULL)
  abline(h = 0, col = "grey")
  return(TRUE)
}
plot_unbiasedness = function(result_, true_, ylim, hypo, effect, confounder, calibration){
  cex.lab = 2
  cex.main = 2.5
  cex.axis = 1.5

  if(confounder == FALSE){
    main = paste("No confounding \n", effect, ', ', hypo, sep = '')
  }else{
    if(calibration == FALSE){
      main = paste("Confounding, unadjusted \n", effect, ', ', hypo, sep = '')
    }else{
      main = paste("Confounding, adjusted \n", effect, ', ', hypo, sep = '')
    }
  }
  xlab = "Time"
  if(effect == "DE"){
    ylab = expression(Delta[DE](t))
  }else{
    ylab = expression(Delta[IE](t))
  }

  time_axis = seq(0, 2.5, 0.01)
  index_all = NULL
  ave_DE = rep(0, length(time_axis))
  plot(NULL, xlim = c(0, 3), ylim = ylim, xlab = xlab, ylab = '', main = main, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
  title(ylab = ylab, line = 2.2, cex.lab = cex.lab)
  for(i in 1:length(result_)){
    ave_DE = ave_DE + approx(result_[[i]]$time, result_[[i]]$effect, xout = time_axis, method = 'constant', rule = 2)$y
    lines(result_[[i]]$time, result_[[i]]$effect, col = 'grey', lwd = .5, type = 's')
    index_all = c(index_all, floor(result_[[i]]$sick * 1000))
  }
  index_all = sort(index_all)
  time_slot = approx(x = index_all, y = 1:length(index_all), xout = unique(index_all), ties = "max", rule = 2, method = "constant")
  time_slot$y = c(time_slot$y[1], diff(time_slot$y))

  small_line = c(ylim[1], ylim[1] + (ylim[2] - ylim[1]) / 15)
  for(i in 1:length(time_slot$x)){
    lines(rep(0.001 * i, 2), small_line, lwd = time_slot$y[i]/max(time_slot$y))
  }
  legend(x = 0, y = ylim[1] + 4 * (ylim[2] - ylim[1]) / 15, legend = c("Average", "True Value"), cex = 1.5, lty = c(1, 2))

  ave_DE = ave_DE/length(result_)
  lines(time_axis, ave_DE, type = 's', lty = 1, lwd = 2)
  lines(true_$time, true_$hazard, type = 's', lty = 2, lwd = 2)
}
plot_sen_ana = function(result){
  width = 400
  height = 400
  for(i in 1:5){
    now_col = rep(0, 11)
    now_gamma = result$DE$sensitivity_analysis[[i]]$gamma
    now_effect = result$DE$sensitivity_analysis[[i]]$effect
    now_col[1] = getcolor(now_gamma[1])
    lwd = 1 + (now_gamma[1] == 0)

    # png(file = paste("/Users/js/Desktop/CHH2020/DE_", i, ".png", sep = ''), width = width, height = height)
    plot(now_effect[[1]]$time/365.25, now_effect[[1]]$effect, type = 's', ylim = c(-0.05, 0.15), col = now_col[1], xlab = "", ylab = "",
         main = paste("Sensitivity analysis, case", i, "\n direct effect"), lwd = lwd, cex.main = 2)
    title(xlab = "Time (years)", ylab = expression(Delta[DE](t)), line = 2.25, cex.lab = 1.5)
    for(j in 2:11){
      now_col[j] = getcolor(now_gamma[j])
      lwd = 1 + (now_gamma[j] == 0)
      lines(now_effect[[j]]$time/365.25, now_effect[[j]]$effect, type = 's', col = now_col[j], lwd = lwd)
    }
    abline(h = 0)
    legend("topleft", legend = now_gamma, col = now_col, lty = 1, lwd = 2, title = expression(gamma), ncol = 2)
    # dev.off()




    now_col = rep(0, 11)
    now_gamma = result$IE$sensitivity_analysis[[i]]$gamma
    now_effect = result$IE$sensitivity_analysis[[i]]$effect
    now_col[1] = getcolor(now_gamma[1])
    lwd = 1 + (now_gamma[1] == 0)

    # png(file = paste("/Users/js/Desktop/CHH2020/IE_", i, ".png", sep = ''), width = width, height = height)
    plot(now_effect[[1]]$time/365.25, now_effect[[1]]$effect, type = 's', ylim = c(0, 0.2), col = now_col[1], xlab = "", ylab = "",
         main = paste("Sensitivity analysis, case", i, "\n indirect effect"), lwd = lwd, cex.main = 2)
    title(xlab = "Time (years)", ylab = expression(Delta[IE](t)), line = 2.25, cex.lab = 1.5)
    for(j in 2:11){
      now_col[j] = getcolor(now_gamma[j])
      lwd = 1 + (now_gamma[j] == 0)
      lines(now_effect[[j]]$time/365.25, now_effect[[j]]$effect, type = 's', col = now_col[j], lwd = lwd)
    }
    abline(h = 0)
    legend("topleft", legend = now_gamma, col = now_col, lty = 1, lwd = 2, title = expression(gamma), ncol = 2)
    # dev.off()
  }
}
getcolor = function(gamma){
  if(gamma == 0){
    return("#000000")
  }else{
    gamma255 = floor(abs(gamma) * 255)
    gammacolor = toupper(as.hexmode(gamma255))
    if(gamma <= 0){
      return(paste("#0000", gammacolor, sep = ''))
    }else{
      return(paste("#", gammacolor, "0000", sep = ''))
    }
  }
}

# last edit at 2021/03/25 16:22
