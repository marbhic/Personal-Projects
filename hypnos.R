CGMS2DayByDay <- function(data, dt0 = NULL, inter_gap = 45, tz = ""){
    
    data = data[complete.cases(data),]
    
    ns = length(unique(data$id))
    
    if (ns > 1){
        id = NULL
        rm(list = c("id"))
        
        first = unique(data$id)[1]
        data = data %>% dplyr::filter(id == first)
        warning(paste("Data contains more than 1 subject. Only the first subject with id", first,  "is used for output."))
    }
    
    ### Get glycemic data
    g = as.numeric(data$gl) # in case it's not
    
    ### Get time data
    if (lubridate::is.POSIXct(data$time)){ # Check if already in date format
        tr = data$time
    }else{
        tr = as.character(data$time)
        tr = as.POSIXct(tr, format='%Y-%m-%d %H:%M:%S', tz = tz)
        # Check if any NAs from conversion, this happens if wrong time format (e.g. 25:00:00) or wrong time zone which will affect daylight savings time
        if (any(is.na(tr))){
            warning(paste("During time conversion,", sum(is.na(tr)), "values were set to NA. Check the correct time zone specification."))
            g = g[!is.na(tr)]
            tr = tr[!is.na(tr)]
        }
    }
    
    timeindex = 2:length(tr)
    timediff = difftime(tr[timeindex], tr[timeindex - 1], units = "mins")
    
    ### Check for time sorting
    if (min(timediff) < 0){
        warning(paste("The times for subject", unique(data$id), "are not in increasing order! The times will be sorted automatically."))
        index = order(tr)
        tr = tr[index]
        g = g[index]
        timediff = difftime(tr[timeindex], tr[timeindex - 1], units = "mins")
    }
    
    ### Automatically identify grid width dt0
    if (is.null(dt0)){
        dt0 = as.double(round(median(timediff, na.rm = TRUE)))
    }
    
    if (dt0 > inter_gap){
        stop(paste("Identified measurements frequency,", dt0, "min, is above the maximal interpolation gap of", inter_gap, "min!"))
    }
    
    ### Check for misaligned grid length dt0 across days, and adjust if needed
    if (1440 %% dt0 != 0){
        if (dt0 > 20){ # Larger grid lengths are set to 20 min
            dt0 = 20
        }else{ # Smaller grid lengths are rounded so that they are divisible by 5 min
            remainder = dt0 %% 5
            if (remainder > 2){
                dt0 = dt0 + 5 - remainder
            }else{
                dt0 = dt0 - remainder
            }
        }
    }
    
    ### Create ideal grid to interpolate over, from minimal to maximal day
    ndays = ceiling(as.double(difftime(max(tr), min(tr), units = "days")) + 1)
    dti = rep(dt0, ndays * 24 * 60 /dt0)
    dti_cum = cumsum(dti)
    dti_cum = lubridate::minutes(dti_cum)
    
    # Set up starting point at 00:00am on the first day
    minD = min(tr) # 1st day of measurements
    lubridate::hour(minD) = 0 # set to 00am
    lubridate::minute(minD) = 0 # set to 00:00
    lubridate::second(minD) = 0 # set to 00:00:00
    
    # Create a set of time points for interpolation
    time_out = minD + dti_cum
    
    ### Interpolate on the ideal grid
    new <- as.data.frame(stats::approx(x = tr, y = g, xout = time_out))
    
    ### Adjust to that there is no interpolation between values > inter_gap appart
    ### Thanks to weather_interp function from weathercan R package for inspiration
    inter_gap <- lubridate::minutes(inter_gap)
    timediff <- lubridate::minutes(round(timediff))
    which_gap <- tr[c(timediff > inter_gap, FALSE)]
    missing <- lubridate::interval(which_gap + 1, which_gap + timediff[timediff > inter_gap] - 1)
    missing <- vapply(new$x, FUN.VALUE = TRUE, FUN = function(x, missing) {
        any(lubridate::`%within%`(x, missing))
    }, missing = missing)
    new$y[missing] <- NA
    
    # Next, from ti remove all the ones that are more than dt0 min away from t0
    gd2d = matrix(new$y, nrow = ndays, byrow = TRUE)
    
    # Assign rownames that correspond to actual dates
    actual_dates = as.Date(minD) + lubridate::days(0:(ndays - 1))
    
    return(list(gd2d = gd2d, actual_dates = actual_dates, dt0 = dt0))
}
all_metrics <- function(data){
    # Mean, Median, and Quantile Metrics not included. Summary covers all
    out = list("ADRR" = adrr(data),
               "AUC" = auc(data),
               "CONGA" = conga(data),
               "CV_GLU" = cv_glu(data),
               "CV_Measures" = cv_measures(data),
               "eA1C" = ea1c(data),
               "GMI" = gmi(data),
               "GRADE" = grade(data),
               "GRADE_Euglycemia" = grade_eugly(data),
               "GRADE_Hyperglycemia" = grade_hyper(data),
               "GRADE_Hypoglycemia" = grade_hypo(data),
               "GVP" = gvp(data),
               "HBGI" = hbgi(data),
               "LBGI" = lbgi(data),
               "Hyper_Index" = hyper_index(data),
               "Hypo_Index" = hypo_index(data),
               "IGC" = igc(data),
               "IQR_GLU" = iqr_glu(data),
               "J_Index" = j_index(data),
               "M_Value" = m_value(data),
               "Mad_GLU" = mad_glu(data),
               "MAGE" = mage(data),
               "MODD" = modd(data),
               "Percent_Above" = above_percent(data),
               "Percent_Below" = below_percent(data),
               "Percent_In_Range" = in_range_percent(data),
               "Range" = range_glu(data),
               "SD_GLU" = sd_glu(data),
               "SD_Measures" = sd_measures(data),
               "SD_ROC" = sd_roc(data),
               "Summary" = summary_glu(data))
    outTable <- out %>%
        Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="id"), .)
    
    return(outTable)
}
 

dexcom1 <- as.data.frame(split_data["dexcom"])
org <- split(dexcom, dexcom$id=="70497")
org1 <- as.data.frame(org["TRUE"])
View(org1)
err1 <-split(dexcom,dexcom$id=="70616")
err1 <- as.data.frame(err1["TRUE"])
err2 <-split(abbott,abbott$id=="70417")
err2 <-as.data.frame(err2["TRUE"])

#Separates cgm data by device
split_data <- split(cgm, cgm$device)
dexcom <- as.data.frame(split_data["dexcom"])
abbott <- as.data.frame(split_data["abbott"])

#changes column names so it can be recognized by iglu functions
names(dexcom)[names(dexcom) == "dexcom.id"] <- "id"
names(dexcom)[names(dexcom) == "dexcom.time"] <- "time"
names(dexcom)[names(dexcom) == "dexcom.gl"] <- "gl"
names(abbott)[names(abbott) == "abbott.time"] <- "time"
names(abbott)[names(abbott) == "abbott.gl"] <- "gl"
names(abbott)[names(abbott) == "abbott.id"] <- "id"

#Runs all_metrics function
abbott_cgm <- all_metrics(abbott)
dexcom_cgm <- all_metrics(dexcom)

#Changes all_metrics output into numeric dataframe
abbott_cgm <- as.data.frame(abbott_cgm)
abbott_cgm$id <- as.numeric(abbott_cgm1$id)
dexcom_cgm <- as.data.frame(dexcom_cgm)
dexcom_cgm$id <- as.numeric(dexcom_cgm$id)

#Creates correlation matrix
cor_abbott <- cor(abbott_cgm, use= "complete.obs",
                  method = c("pearson", "kendall", "spearman"))
cor_dexcom <- cor(dexcom_cgm, use= "complete.obs",
                  method = c("pearson", "kendall", "spearman"))

colors <- c("#67001F", "#B2182B", "#D6604D", "#F4A582",
            "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
            "#4393C3", "#2166AC", "#053061")
colors <- colors[length(colors):1]
col2 <- colorRampPalette(colors)

#creates correlation plots
par(mfrow = c(1,1))
corrplot(cor_dexcom, col = col2(20), method = "color", 
         type = "lower",tl.srt = 40, title= "Dexcom", 
         diag = F, tl.cex = 0.6, cl.cex = 0.6, mar=c(0,0,2,0))
corrplot(cor_abbott, col = col2(20), method = "color", 
         type = "lower",tl.srt = 40, title= "Abbott", 
         diag = F, tl.cex = 0.6, cl.cex = 0.6, mar=c(0,0,2,0))
dev.off()

#registry data
registry <- read_excel("registry.xls")
registry_data <- as.data.frame(registry[1:3])
registry_data$BMI = registry$BMI

#add registry data to cgm data
abbott_reg <- merge(abbott_cgm,registry_data, by = "id")
dexcom_reg <- merge(dexcom_cgm,registry_data, by = "id")

#correlation plot for merged data
cor_abbott_reg <- cor(abbott_reg, use= "complete.obs",
                  method = c("pearson", "kendall", "spearman"))
cor_dexcom_reg <- cor(dexcom_reg, use= "complete.obs",
                  method = c("pearson", "kendall", "spearman"))
corrplot(cor_dexcom_reg, col = col2(20), method = "color", 
         type = "lower",tl.srt = 40, title= "Dexcom", 
         diag = F, tl.cex = 0.6, cl.cex = 0.6, mar=c(0,0,2,0))
corrplot(cor_abbott_reg, col = col2(20), method = "color", 
         type = "lower",tl.srt = 40, title= "Abbott", 
         diag = F, tl.cex = 0.6, cl.cex = 0.6, mar=c(0,0,2,0))
dev.off()


# Do heatmap
abbott_mat = as.matrix(abbott_cgm[ , -1])
rownames(abbott_mat) = as.character(abbott$id)

# Do centering and scaling of all metrics before drawing the heatmap
abbott_mat_scale = scale(abbott_mat)

# Rename the metrics to make the plots nices
metric_names = colnames(abbott_mat_scale)
metric_names[metric_names == "adrr"]="ADRR"
metric_names[metric_names == "hourly_auc"]="AUC"
metric_names[metric_names == "conga"]="CONGA"
metric_names[metric_names == "cv"]="CV"
metric_names[metric_names == "ea1c"]="eA1C"
metric_names[metric_names == "gmi"]="GMI"
metric_names[metric_names == "grade_eugly"]="GRADE eugly"
metric_names[metric_names == "grade_hyper"]="GRADE hyper"
metric_names[metric_names == "grade_hypo"]="GRADE hypo"
metric_names[metric_names == "grade"]="GRADE"
metric_names[metric_names == "gvp"]="GVP"
metric_names[metric_names == "hbgi"]="HBGI"
metric_names[metric_names == "lbgi"]="LBGI"
metric_names[metric_names == "hyper_index"]="Hyper Index"
metric_names[metric_names == "hypo_index"]="Hypo Index"
metric_names[metric_names == "above_140"]="% above 140"
metric_names[metric_names == "above_180"]="% above 180"
metric_names[metric_names == "above_200"]="% above 200"
metric_names[metric_names == "above_250"]="% above 250"
metric_names[metric_names == "below_50"]="% below 50"
metric_names[metric_names == "below_80"]="% below 80"
metric_names[metric_names == "range"]="Range"
metric_names[metric_names == "iqr"]="IQR"
metric_names[metric_names == "igc"]="IGC"
metric_names[metric_names == "j_index"]="J index"
metric_names[metric_names == "m_value"]="M value"
metric_names[metric_names == "mage"]="Mage"
metric_names[metric_names == "1st Qu."]="1st quartile"
metric_names[metric_names == "3rd Qu."]="3rd quartile"
metric_names[metric_names == "in_range_70_140"]="% in 70-140"
metric_names[metric_names == "in_range_70_180"]="% in 70-180"
metric_names[metric_names == "in_range_80_200"]="% in 80-200"
metric_names[metric_names == "modd"]="MODD"
metric_names[metric_names == "sd"]="Sd"
metric_names[metric_names == "sd_roc"]="Sd ROC"
metric_names[metric_names == "Min."]="Min"
metric_names[metric_names == "Max."]="Max"
colnames(abbott_mat_scale) = metric_names
p = pheatmap(t(abbott_mat_scale), cutree_rows = 6, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
             clustering_method = "complete", angle_col = 0)
p
