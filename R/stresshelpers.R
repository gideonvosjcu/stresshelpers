# Investigating Wearable Sensor Biomarkers for Chronic Stress Measurement and Analysis
# Gideon Vos, Master of Philosophy, James Cook University
# 14 September 2021

# Citations:

# WESAD (Wearable Stress and Affect Detection)
# Philip Schmidt, Attila Reiss, Robert Duerichen, Claus Marberger, and Kristof Van Laerhoven. 2018.
# Introducing WESAD, a Multimodal Dataset for Wearable Stress and Affect Detection.
# In Proceedings of the 20th ACM International Conference on Multimodal Interaction (ICMI '18).
# Association for Computing Machinery, New York, NY, USA, 400–408. DOI:https://doi.org/10.1145/3242969.3242985

# AffectiveROAD:
# Neska El Haouij, Jean-Michel Poggi, Sylvie Sevestre-Ghalila, Raja Ghozi, and Mériem Jaïdane.
# 2018. AffectiveROAD system and database to assess driver's attention. In Proceedings of the 33rd Annual ACM Symposium on Applied Computing (SAC '18).
# ACM, New York, NY, USA, 800-803. DOI: https://doi.org/10.1145/3167132.3167395. https://dam-prod.media.mit.edu/x/2021/06/14/AffectiveROAD_Data_w1dqSB9.zip

# The SWELL Knowledge Work Dataset for Stress and User Modeling Research
# Koldijk, S., Sappelli, M., Verberne, S., Neerincx, M., & Kraaij, W. (2014).
# The SWELL Knowledge Work Dataset for Stress and User Modeling Research.
# To appear in: Proceedings of the 16th ACM International Conference on Multimodal Interaction (ICMI 2014) (Istanbul, Turkey, 12-16 November 2014).
# The dataset can be accessed medio 2015 here: http://persistent-identifier.nl/?identifier=urn:nbn:nl:ui:13-kwrv-3e.

# Non-EEG Dataset for Assessment of Neurological Status
# Birjandtalab, Javad, Diana Cogan, Maziyar Baran Pouyan, and Mehrdad Nourani,
# A Non-EEG Biosignals Dataset for Assessment and Visualization of Neurological Status,
# 2016 IEEE International Workshop on Signal Processing Systems (SiPS), Dallas, TX, 2016, pp. 110-114. doi: 10.1109/SiPS.2016.27

# Toadstool: A Dataset for Training Emotional Intelligent Machines Playing Super Mario Bros
# Svoren, H., Thambawita, V., Halvorsen, P., Jakobsen, P., Garcia-Ceja, E., Noori, F. M., … Hicks, S. (2020, February 28).
# https://doi.org/10.31219/osf.io/4v9mp

# K-EmoCon, a multimodal sensor dataset for continuous emotion recognition in naturalistic conversations
# Cheul Young Park, Narae Cha, Soowon Kang, Auk Kim, Ahsan Habib Khandoker, Leontios Hadjileontiadis, Alice Oh, Yong Jeong, & Uichin Lee. (2020).
# In Scientific Data (1.0.0, Vol. 7, Number 1, p. 293). Zenodo. https://doi.org/10.5281/zenodo.3931963

# Multilevel Monitoring of Activity and Sleep in Healthy People (version 1.0.0)
# Rossi, A., Da Pozzo, E., Menicagli, D., Tremolanti, C., Priami, C., Sirbu, A., Clifton, D., Martini, C., & Morelli, D. (2020).. PhysioNet.
# https://doi.org/10.13026/cerq-fc86.

#########################################################################################################################################################
# Helper functions
#########################################################################################################################################################

# function to compute total within-cluster sum of square

#' Compute Total Within-Cluster Sum of Square
#'
#' @param data The data frame to apply kmeans to
#' @param k The starting number of clusters
#' @return The total within-cluster sum of sqaure
#' @export
wss <- function(data, k) {
  kmeans(data, k, nstart = 5 )$tot.withinss
}

#' Predict clusters on new data given existing kmeans model
#'
#' @param object The kmeans model
#' @param newdata The new data frame
#' @param method Centers or Classes
#' @return A vector of predictions
#' @export
predict.kmeans <- function(object, newdata, method = c("centers", "classes")) {
  method <- match.arg(method)
  centers <- object$centers
  ss_by_center <- apply(centers, 1, function(x) {
    colSums((t(newdata) - x) ^ 2)
  })
  best_clusters <- apply(ss_by_center, 1, which.min)
  if (method == "centers") {
    centers[best_clusters, ]
  } else {
    best_clusters
  }
}

#' Get the mode of a vector
#'
#' @param v The vector
#' @return The mode
#' @export
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#########################################################################################################################################################
# Cortisol modelling functions
#########################################################################################################################################################

#' Connects EDA peaks
#'
#' @param peaks The EDA vector
#' @param peak_height EDA cut-off height to consider as a peak
#' @param window_size Window size of EDA peak
#' @return The adjusted vector
#' @export
connect_peaks <- function(peaks, peak_height, window_size)
{
  new_peaks <- peaks
  series_length <- length(peaks)
  for (index in (1:length(peaks)))
  {
    subset <- peaks[index:(index+window_size)]
    subset[is.na(subset)] <- 0
    if (max(subset) >= peak_height)
    {
      subset[which(subset >= peak_height)] <- max(subset)
    }
    new_peaks[index:(index+window_size)] <- subset
  }
  return (new_peaks[1:series_length])
}

#' Model a Cortisol curve
#'
#' @param x The EDA vector
#' @param peak_height EDA cut-off height to consider as a peak
#' @return A vector of Cortisol measurements
#' @export
model_cortisol <- function(x, peak_height)
{
  peak_last_time <- 90 * 60 # 90 minutes from EDA peak
  cortisol <- rep(2.5, length(x))
  max_eda <- max(x) # determine scale

  index <- 1
  while (index < length(x))
  {
    if (x[index] > peak_height)
    {
      subset <- x[index:(index+peak_last_time-1)]
      # from index to peak, scale up
      max_peak_height <- 20 /max_eda * x[index]
      step_size <- (max_peak_height - 2.5) / (15 * 60) # step up from 2.5 to max for 15 minutes
      step <- step_size + 2.5
      for (j in 1:(15*60))
      {
        cortisol[index+j-1] <- step
        step <- step + step_size
      }
      step_size <- (max_peak_height - 2.5) / (75 * 60) # step up from 2.5 to max for 15 minutes
      step <- step_size
      for (j in 1:(75 * 60))
      {
        cortisol[index+j-1+ (15 * 60)] <- max_peak_height - step
        step <- step + step_size
      }

      index <- index + 25 # 25 second window
    }
    else
    {
      index <- index + 1
    }
  }
  cortisol <- cortisol[1:length(x)]
  return (cortisol)
}
#########################################################################################################################################################
# Empatica E4 data reader - courtesy of https://github.com/bwrc/empatica-r
#########################################################################################################################################################
organise_data <- function(data, samplingrate = NULL) {
  out <- vector(mode = "list", length = ncol(data))
  if (!is.null(samplingrate))
  {
    time_vector <- seq.int(0, (nrow(data) - 1)) / samplingrate
  }
  else
  {
    time_vector <- NULL
  }
  for (i in seq.int(ncol(data))) {
    signal_name <- colnames(data)[i]
    out[[i]]$data <- data[,i]
    out[[i]]$t <- time_vector
    out[[i]]$samplingrate <- samplingrate
  }
  names(out) <- colnames(data)
  return (out)
}

read_empatica_ibi <- function(filename, signal_names = c("time", "ibi")) {
  t_start <- scan(filename, nmax = 1, what = numeric(), nlines = 1, skip = 0, sep = ",", quiet = TRUE)
  data <- read.csv(filename, header = FALSE, sep = ",", skip = 1)
  if (!is.null(signal_names))
  {
    colnames(data) <- signal_names
  }
  out <- organise_data(data)
  out[["ibi"]][["t"]] <- out[["time"]][["data"]]
  out[["ibi"]][["samplingrate"]] <- NULL
  out[["time"]] <- NULL
  return(out)
}

read_header <- function(filename) {
  scan(filename, nmax = 1, what = numeric(), nlines = 1, skip = 0, sep = ",", quiet = TRUE)
}

read_empatica_events <- function(f) {
  tmp <- scan(f, quiet = TRUE)
  if (length(tmp) > 0)
  {
    return (data.frame("id" = seq.int(length(tmp)), "time_raw"  = tmp,
                       "timestamp" = as.POSIXct(tmp, tz = "GMT", origin = "1970-01-01"), "timedelta" = 0))
  }
  else
  {
    return (data.frame("id" = numeric(), "time_raw" = numeric(), "timestamp" = as.POSIXct(character()), "timedelta" = numeric()))
  }
}

read_empatica_gen <- function(filename, signal_names = NULL) {
  t_start <- scan(filename, nmax = 1, what = numeric(), nlines = 1, skip = 0, sep = ",", quiet = TRUE)
  samplingrate <- scan(filename, nmax = 1, what = numeric(), nlines = 1, skip = 1, sep = ",", quiet = TRUE)
  data <- read.csv(filename, header = FALSE, sep = ",", skip = 2)

  if (!is.null(signal_names))
  {
    colnames(data) <- signal_names
  }
  return (organise_data(data, samplingrate))
}

new_recording <- function() {
  recording <- list()
  recording$properties <- list()
  recording$signal <- list()
  recording$events <- list()
  recording$properties$header <- list()
  recording$properties$time.start.raw <- NA
  recording$properties$time.start <- NA
  recording$properties$time.stop.raw <- NA
  recording$properties$time.stop <- NA
  recording$properties$subject <- NA
  recording$properties$format <- NA
  recording$properties$format.long <- NA
  recording$properties$device.type <- NA
  recording$properties$device.version <- NA
  recording$properties$length <- NA
  return (recording)
}

#' Reads an Empatica E4 folder of data into a data frame
#'
#' @param exdir The folder containing the E4 data
#' @param header.only Header only, T or F
#' @return A data frame containing the data
#' @export
read.empatica <- function(exdir, header.only = FALSE) {
  recording <- new_recording()
  filelist <- list.files(exdir, pattern = "*.csv", full.names = TRUE)

  out <- list()
  for (f in filelist) {
    signal_name <- tolower(gsub(".csv", "", basename(f)))
    signal_names <- NULL
    if (signal_name == "acc")
      signal_names <- c("acc_x", "acc_y", "acc_z")
    if (signal_name == "bvp")
      signal_names <- c("bvp")
    if (signal_name == "eda")
      signal_names <- c("eda")
    if (signal_name == "hr")
      signal_names <- c("hr")
    if (signal_name == "temp")
      signal_names <- c("temp")
    if (signal_name == "ibi") {
      signal_names <- c("ibi")
      out <- c(out, read_empatica_ibi(f))
      signal_names <- NULL
    }
    if (signal_name == "tags") {
      recording$events <- read_empatica_events(f)
      signal_names <- NULL
    }

    if (! is.null(signal_names))
      out <- c(out, read_empatica_gen(f, signal_names))
  }
  header <- read_header(filelist[[1]])
  recording$signal <- out
  if (length(names(recording$signal)) > 1)
  {
    recording$properties$length <- rev(recording$signal[[names(recording$signal)[1]]]$t)[1]
  }
  else
  {
    recording$properties$length <- 0
  }
  recording$properties$time.start.raw <- header
  recording$properties$time.start <- as.POSIXct(header, tz = "GMT", origin = "1970-01-01")
  recording$properties$time.stop.raw <- header + recording$properties$length
  recording$properties$time.stop <- as.POSIXct(header, tz = "GMT", origin = "1970-01-01") + as.difftime(recording$properties$length, units = "secs")
  if (length(recording$events) > 0)
    recording$events$timedelta <- as.numeric(difftime(recording$events$timestamp, recording$properties$time.start, units = "secs"))
  recording$properties$device.type <- "Empatica"
  recording$properties$device.version <- "E4"
  return (recording)
}

#########################################################################################################################################################
# PhysioNet file reader - courtesy of https://rdrr.io/github/Absox/wfdb.R/f/
#########################################################################################################################################################
wfdbdesc <- function(header.filename, read.datetime = TRUE) {
  # Get root directory of header file, as we may need to read other files using that path
  root.dir <- ""
  if (grepl("/", header.filename)) {
    root.dir <- substr(header.filename, 1, regexpr("/[^/]*$", header.filename))
    header.filename <- substr(header.filename, regexpr("/[^/]*$", header.filename)+1, nchar(header.filename))
  }
  # If filename lacks extension, add extension
  if (!grepl(".hea$", header.filename)) {
    header.filename <- sprintf("%s.hea", header.filename)
  }

  # Read in header file for parsing
  header.filepath <- sprintf("%s%s", root.dir, header.filename)
  header.file <- file(header.filepath)
  header.data <- trimws(strsplit(readChar(header.file,file.size(header.filepath)),"\n")[[1]])
  close(header.file)

  # Process the first line of the header (record line):
  record.line <- strsplit(header.data[1]," ")[[1]]
  # Determine what type of record this is: single segment or multi-segment, number of signals, sampling frequency,
  # and record length
  record.name <- record.line[1]
  num.signals <- as.numeric(record.line[2])
  if (grepl("/", record.line[3])) {
    sample.frequency <- as.numeric(substr(record.line[3],1,regexpr("/",record.line[3])-1))
  } else {
    sample.frequency <- as.numeric(record.line[3])
  }
  record.length <- as.numeric(record.line[4])

  # Single segment datetimes
  if (read.datetime) {
    base.datetime <- lubridate::dmy(record.line[6]) + lubridate::hms(record.line[5])
  } else {
    base.datetime <- NA
  }

  if (grepl("/", record.name)) {
    # Multi segment file
    record.type <- "multi"
    num.segments <- as.numeric(substr(record.name, regexpr("/",record.name)+1,nchar(record.name)))
    record.name <- substr(record.name, 1, regexpr("/",record.name)-1)
    # Read segment names and segment lengths
    segments <- strsplit(header.data[2:(num.segments+1)], " ")
    segment.lengths <- sapply(segments, function(x) as.numeric(x[2]))
    segment.names <- sapply(segments, function(x) x[1])

    # From first valid segment, read descriptions, adc gains, and adc zeros
    first.segment <- segment.names[min(which(segment.names != "~"))]
    layout.filepath <- sprintf("%s%s.hea", root.dir, first.segment)
    layout.file <- file(layout.filepath)
    layout.data <- trimws(strsplit(readChar(layout.file, file.size(layout.filepath)), "\n")[[1]])
    close(layout.file)

    signal.lines <- strsplit(layout.data[2:(num.signals+1)], " ") # TODO
    signal.format <- NULL
  } else {
    # Single segment, can proceed to read signal specification from current file
    record.type <- "single"
    signal.lines <- strsplit(header.data[2:(num.signals+1)], " ")
    segment.names <- signal.lines[[1]][1];
    segment.lengths <- record.length
    signal.format <- as.numeric(signal.lines[[1]][2])
  }

  # Read signal descriptions, adc gains, units of measure, and zeros from signal specification
  descriptions <- sapply(signal.lines, function(x) {
    if (length(x) > 9) {
      return(paste(x[-(1:8)],collapse=" "))
    } else {
      return(x[9])
    }
  })
  # This extra code determines if there is a baseline
  adc.gains <- sapply(signal.lines, function(x) {
    if (grepl("\\(.*\\)", x[3])) {
      return(as.numeric(substr(x[3], 1, regexpr("\\(", x[3])-1)))
    } else {
      return(as.numeric(substr(x[3], 1, regexpr("/", x[3])-1)))
    }

  })
  adc.units <- sapply(signal.lines, function(x) substr(x[3], regexpr("/", x[3])+1, nchar(x[3])))
  # If no baseline is present, zero is set to ADC zero. Otherwise it is equal to the baseline.
  adc.zeros <- sapply(signal.lines, function(x) {
    if (grepl("\\(.*\\)", x[3])) {
      return(as.numeric(substr(x[3], regexpr("\\(", x[3])+1, regexpr("\\)", x[3])-1)))
    } else {
      return(as.numeric(x[5]))
    }
  })

  return(list(root.dir = root.dir, record.name = record.name, record.type = record.type,
              sample.frequency = sample.frequency,  record.length = record.length, base.datetime = base.datetime,
              signal.format = signal.format, num.signals = num.signals, descriptions = descriptions, adc.gains = adc.gains,
              adc.zeros = adc.zeros, adc.units = adc.units, segment.names = segment.names,
              segment.lengths = segment.lengths))
}

#' Reads a a PhysioNET file
#'
#' @param header.filename The header file
#' @return A data frame containing the data
#' @export
rdsamp <- function(header.filename) {
  siginfo <- wfdbdesc(header.filename, FALSE)

  read.data <- function(data.filename, signal.format, num.signals, signal.length, adc.gains, adc.zeros) {
    data.file = file(data.filename, "rb")

    if (signal.format == 16) {
      signal.data <- t((array(readBin(data.file, integer(), n = file.size(data.filename)/2, size = 2,
                                      endian = "little"), dim=c(num.signals, signal.length)) - adc.zeros)/adc.gains)
    } else if (signal.format == 80) {
      raw.data <- readBin(data.file, integer(), n = file.size(data.filename), size = 1, endian = "little", signed = FALSE)
      signal.data <- t((array(raw.data - 128, dim=c(num.signals, signal.length)) - adc.zeros)/adc.gains)
    }

    close(data.file)
    return(signal.data)
  }

  if (siginfo$record.type == "single") {
    # If single-segment file, we can simply read and return the first segment.
    return (read.data(data.filename = sprintf("%s%s", siginfo$root.dir, siginfo$segment.names),
                      signal.format = siginfo$signal.format, num.signals = siginfo$num.signals,
                      signal.length = siginfo$record.length, adc.gains = siginfo$adc.gains, adc.zeros = siginfo$adc.zeros))
  } else {
    # We have to allocate a table for the unified data, read each segment, and place it into its proper place
    signal.data <- array(NA, dim = c(siginfo$record.length, siginfo$num.signals))
    segment.starts <- c(cumsum(c(1,siginfo$segment.lengths[1:(length(siginfo$segment.lengths)-1)])))
    segment.ends <- c(segment.starts[-1]-1, siginfo$record.length)
    segment.starts <- segment.starts[siginfo$segment.names != "~" & siginfo$segment.lengths > 0]
    segment.ends <- segment.ends[siginfo$segment.names != "~" & siginfo$segment.lengths > 0]
    # Recursively read data from segments
    valid.segments = sapply(siginfo$segment.names[siginfo$segment.names != "~" & siginfo$segment.lengths > 0],
                            function(x) sprintf("%s%s", siginfo$root.dir, x), USE.NAMES = FALSE)

    # Need to check that there are any valid segments
    if(length(valid.segments) > 0) {
      segment.headers <- lapply(valid.segments, function(x) wfdbdesc(x, FALSE))
      segment.data <- lapply(valid.segments, rdsamp)
      # Yo dawg, I heard you like lambdas, so I put a lambda in my lambda so you can lambda while you lambda
      indices <- lapply(segment.headers, function(x) sapply(x$descriptions,
                                                            function (y) which(siginfo$descriptions == y)))
      # Unify the data
      for (c in 1:length(segment.data)) {
        signal.data[segment.starts[c]:segment.ends[c],indices[[c]]] = segment.data[[c]]
      }
    }

    return(signal.data)
  }
}

#########################################################################################################################################################
# Feature engineering routines
#########################################################################################################################################################

#' Downsample a series to match size of another series given by rate
#'
#' @param data A vector of sensor data
#' @param rate The sampling rate to downsize to
#' @return A downsampled vector
#' @export
downsample <- function(data, rate)
{
  result <- c()
  for(i in seq(1, length(data), by = rate))
  {
    result <- c(result, mean(data[i:(i+rate)]))
  }
  return (result)
}

# create windowed features of a specific window size
# does not slide

#' Creates new features from Empatica E4 data
#'
#' @param data A data frame of EDA and HR sensor data
#' @param window_size Window size to use for rolling windows
#' @return An expanded data frame of features
#' @export
rolling_features <- function(data, window_size)
{
  for (row in seq(1,nrow(data), by=window_size))
  {
    subset <- data[row:(row+(window_size-1)),]
    data[row:(row+(window_size-1)), "edamean"] <- mean(subset$eda)
    data[row:(row+(window_size-1)), "edamedian"] <- median(subset$eda)
    data[row:(row+(window_size-1)), "edastd"] <- sd(subset$eda)
    data[row:(row+(window_size-1)), "edavar"] <- var(subset$eda)
    data[row:(row+(window_size-1)), "edamin"] <- min(subset$eda)
    data[row:(row+(window_size-1)), "edamax"] <- max(subset$eda)
    data[row:(row+(window_size-1)), "edaskew"] <- skewness(subset$eda)
    data[row:(row+(window_size-1)), "edakurt"] <- kurtosis(subset$eda)
    data[row:(row+(window_size-1)), "edarange"] <- max(subset$eda) - min(subset$eda)
    data[row:(row+(window_size-1)), "hrmean"] <- mean(subset$hr)
    data[row:(row+(window_size-1)), "hrmedian"] <- median(subset$hr)
    data[row:(row+(window_size-1)), "hrstd"] <- sd(subset$hr)
    data[row:(row+(window_size-1)), "hrmin"] <- min(subset$hr)
    data[row:(row+(window_size-1)), "hrmax"] <- max(subset$hr)
    data[row:(row+(window_size-1)), "hrvar"] <- var(subset$hr)
    data[row:(row+(window_size-1)), "hrskew"] <- skewness(subset$hr)
    data[row:(row+(window_size-1)), "hrkurt"] <- kurtosis(subset$hr)
    data[row:(row+(window_size-1)), "hrrange"] <- max(subset$hr) - min(subset$hr)
    data[row:(row+(window_size-1)), "cov1"] <- cov(subset$hr, subset$eda)
  }
  data[is.nan(data$edakurt),"edakurt"] <-0
  data[is.nan(data$edaskew),"edaskew"] <-0
  data <- na.omit(data)
  # move metric to be last column
  metric <- data$metric
  data$metric <- NULL
  data$metric <- metric
  return (data)
}
