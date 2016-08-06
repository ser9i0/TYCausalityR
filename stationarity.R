# Test a time series stationarity usin ADF and KPSS test
# Return max order of integration
stationary_test <- function(dataset) 
{
  max_order_integration<-0
  max_order_integration_adf<-0
  diff_order<-0

  # Augmented Dickey-Fuller unit-root test
  # The null hypothesis is non-stationarity
  adf_htest <-adf.test(dataset)
  adf_p <- adf_htest$p.value
  adf_s <- adf_htest$statistic

  if( adf_p < 0.1 )
  {
    # Reject of null hypothesis (stationary) with at least 10% significance level
    message("....ADF Test: is STATIONARY with p-value = ",adf_p," and statistic = ",adf_s)
  }
  else
  {
    # Differencing the time series until stationarity is found
    message("....ADF Test: is NON-STATIONARY with p-value = ",adf_p," and statistic = ",adf_s)
    message("------Calculating order of integration...")
    while( adf_p >= 0.1 )
    {
      diff_order<-diff_order+1
      adf_htest <-adf.test(diff(dataset,1,diff_order))
      adf_p <- adf_htest$p.value
      adf_s <- adf_htest$statistic
    }
    message("------ADF Test: is STATIONARY with p-value = ",adf_p," and statistic = ",adf_s," and order of integration I(",diff_order,")")
    
    if( diff_order > max_order_integration_adf )
      max_order_integration_adf<-diff_order
  }

  max_order_integration_kpss<-0    
  diff_order<-0

  # KPSS (Kwiatkowski–Phillips–Schmidt–Shin) test
  # The null hypothesis is stationarity
  kpss_htest <-kpss.test(dataset)
  kpss_p <- kpss_htest$p.value
  kpss_s <- kpss_htest$statistic

  if( kpss_p < 0.1 )
  {
    message("....KPSS Test: is NON-STATIONARY with p-value = ", kpss_p, " and statistic = ", kpss_s)
    message("------Calculating order of integration...")
    while( kpss_p < 0.1 )
    {
      diff_order<-diff_order+1
      kpss_htest <-kpss.test(diff(dataset,1,diff_order))
      kpss_p <- kpss_htest$p.value
      kpss_s <- kpss_htest$statistic
    }
    message("------KPSS Test: is STATIONARY with p-value = ",kpss_p," and statistic = ",kpss_s," and order of integration I(",diff_order,")")

    if( diff_order > max_order_integration_kpss )
         max_order_integration_kpss<-diff_order
  }
  else
  {
    # Could not reject the null hypothesis (non-stationary) with at least 10% significance level
    message("....KPSS Test: is STATIONARY with p-value = ", kpss_p, " and statistic = ", kpss_s)
  }

  # Given a group of time series, the maximum order of integration (m) is the highest order of integration
  # of each time series. For example if one is I(0) and the other is I(2), then m = 2, etc.

  if( max_order_integration_adf > max_order_integration_kpss )
  {
    max_order_integration<-max_order_integration_adf
    message("WARNING!: Max order of integration differs from ADF Test (",max_order_integration_adf,") to KPSS Test (",max_order_integration_kpss,"). The highest one will be used.")
  }
  else if( max_order_integration_kpss > max_order_integration_adf )
  {
    max_order_integration<-max_order_integration_kpss
    message("WARNING!: Max order of integration differs from ADF Test (",max_order_integration_adf,") to KPSS Test (",max_order_integration_kpss,"). The highest one will be used.")
  }
  else
  {
    max_order_integration<-max_order_integration_adf
    message("....Max order of integration: ",max_order_integration)
  }

  result <- list(max_order_integration = max_order_integration)
  return(result)
}
