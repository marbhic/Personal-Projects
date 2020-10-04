
check_data_columns =  function(data, id = 'id', time = 'time', gl = 'gl', time_check = FALSE, tz = ""){
    if (is.vector(data)) {
        output = as.double(data)
        output = data.frame(gl = output,
                            id = 1,
                            time = NA_real_)
        attr(output, "is_vector") = TRUE
    } else {
        cols_in = c(id, time, gl) %in% names(data)
        if (!all(cols_in)) {
            s = setdiff(c(id, time, gl), names(data))
            msg = paste0("Columns: ", paste0(s, collapse = ", "),
                         " are not present in the data")
            stop(msg)
        }
        if (time_check) {
            if (!lubridate::is.POSIXct(data$time)){ # Check if already in date format
                tr = as.character(data$time)
                data$time = as.POSIXct(tr, format='%Y-%m-%d %H:%M:%S', tz = tz)
            }
        }
        output = data[, c(id, time, gl), drop = FALSE]
        colnames(output) = c("id", "time", "gl")
        attr(output, "is_vector") = FALSE
    }
    return(output)
}
sd_measures1 <- function(data, dt0 = NULL, inter_gap = 45, tz = ""){
    
    id = NULL
    rm(list = c("id"))
    data = check_data_columns(data, time_check = TRUE)
    subject = unique(data$id)
    ns = length(subject)
    
    # Calculate uniform grid for all subjects
    gdall = list()
    for (i in 1:ns){
        if (i != 1){
            dt0 = out$dt0
        }
        out = data %>%
            dplyr::filter(id == subject[i]) %>%
            CGMS2DayByDay(tz = tz, dt0 = dt0, inter_gap = inter_gap)
        gdall[[i]] <- out$gd2d
    }
    dt0 = out$dt0
    
    results = lapply(
        gdall,
        function(gd2d) {
            # SdW - vertical within days
            out = tibble::tibble(id = NA, SdW = mean(apply(gd2d, 1, sd, na.rm = TRUE), na.rm = TRUE))
            # SdHHMM - between time points
            out$SdHHMM = sd(apply(gd2d, 2, mean, na.rm = TRUE), na.rm = TRUE)
            # SdWSH - Within series - for 1 hour window
            win = round(60/dt0) # how many measurements are within 1 hour
            gs = as.vector(t(gd2d))
            #N = length(gs) # total # of measurements
            #ind = rep(1:ceiling(N/win), each = win)[1:N] # consecutive indexing
            #out$SdWSH = mean(tapply(gs, ind, sd, na.rm = TRUE), na.rm = TRUE)
            out$SdWSH = mean(caTools::runsd(gs, k = win, endrule = "trim"), na.rm = TRUE)
            
            # SdDM - "Horizontal" sd
            meanR = apply(gd2d, 1, mean, na.rm = TRUE)
            out$SdDM = sd(meanR, na.rm = TRUE)
            
            # SdB - between days, within timepoints
            out$SdB = mean(apply(gd2d, 2, sd, na.rm = TRUE), na.rm = TRUE)
            
            # SdBDM - between days, within timepoints, corrected for changes in daily means
            med = matrix(rep(meanR, each = ncol(gd2d)), ncol = ncol(gd2d), byrow = TRUE)
            # out$SdBDM = mean(apply(sqrt((gd2d - med)^2), 1, mean, na.rm = TRUE), na.rm = TRUE)
            out$SdBDM = mean(apply(gd2d - med, 2, sd, na.rm = TRUE), na.rm = TRUE)
            
            out
        })
    
    results = dplyr::bind_rows(results)
    results$id = subject
    
    return(results)
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
               "SD_Measures" = sd_measures1(data),
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

na_dexcom <- dexcom[!(dexcom$id == "70616"),]
na_abbott <- abbott[!(abbott$id == "70497"),]

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
abbott_cgm$id <- as.numeric(abbott_cgm$id)
dexcom_cgm <- as.data.frame(dexcom_cgm)
dexcom_cgm$id <- as.numeric(dexcom_cgm$id)

#Creates correlation matrix
cor_abbott <- cor(abbott_cgm, use= "complete.obs",
                  method = c("pearson", "kendall", "spearman"))
cor_dexcom <- cor(dexcom_cgm, use= "complete.obs",
                  method = c("pearson", "kendall", "spearman"))
#test
cor_abbott <- cor(na_abbott, use= "complete.obs",
                  method = c("pearson", "kendall", "spearman"))
cor_dexcom <- cor(na_dexcom, use= "complete.obs",
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


#test
na_dexcom <- dexcom_cgm[!(dexcom_cgm$id == "70616"),]
na_abbott <- abbott_cgm[!(abbott_cgm$id == "70417"),]

# Do heatmap
mecs_mat = as.matrix(na_abbott[ , -1])
rownames(mecs_mat) = as.character(na_abbott$id)

# Do centering and scaling of all metrics before drawing the heatmap
mecs_mat_scale = scale(mecs_mat)

# Rename the metrics to make the plots nices
metric_names = colnames(mecs_mat_scale)
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
colnames(mecs_mat_scale) = metric_names

p = pheatmap(t(mecs_mat_scale), cutree_rows = 6, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
             clustering_method = "complete", angle_col = 0)
p

print(p)

#dexcom
mecs_mat = as.matrix(na_dexcom[ , -1])
rownames(mecs_mat) = as.character(na_dexcom$id)

# Do centering and scaling of all metrics before drawing the heatmap
mecs_mat_scale = scale(mecs_mat)

# Rename the metrics to make the plots nices
metric_names = colnames(mecs_mat_scale)
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
colnames(mecs_mat_scale) = metric_names

p = pheatmap(t(mecs_mat_scale), cutree_rows = 6, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
             clustering_method = "complete", angle_col = 0)
p

#registry data
registry <- read_excel("registry.xls")
registry_data <- as.data.frame(registry[1:3])
registry_data$BMI = registry$BMI

#add registry data to cgm data
abbott_reg <- merge(na_abbott,registry_data, by = "id")
dexcom_reg <- merge(na_dexcom,registry_data, by = "id")

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

#dexcom
mecs_mat = as.matrix(abbott_reg[ , -1])
rownames(mecs_mat) = as.character(abbott_reg$id)

# Do centering and scaling of all metrics before drawing the heatmap
mecs_mat_scale = scale(mecs_mat)

# Rename the metrics to make the plots nices
metric_names = colnames(mecs_mat_scale)
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
colnames(mecs_mat_scale) = metric_names

p = pheatmap(t(mecs_mat_scale), cutree_rows = 6, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
             clustering_method = "complete", angle_col = 0)
p

#dexcom
mecs_mat = as.matrix(dexcom_reg[ , -1])
rownames(mecs_mat) = as.character(dexcom_reg$id)

# Do centering and scaling of all metrics before drawing the heatmap
mecs_mat_scale = scale(mecs_mat)

# Rename the metrics to make the plots nices
metric_names = colnames(mecs_mat_scale)
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
colnames(mecs_mat_scale) = metric_names

p = pheatmap(t(mecs_mat_scale), cutree_rows = 6, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
             clustering_method = "complete", angle_col = 0)
p
