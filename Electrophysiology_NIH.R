pacman::p_load(tidyverse)
data<-read_csv("./paired_pulse_df.csv")
head(data)

ggplot(data,aes(x=time, y=sweep_1) +
  geom_point() +
  xlim(0.05,0.08) +
  ylim(-0.04,0.05))
  
#sweep_of_interest
#EG FXN
#To do
 #1 Find peaks and valleys
 #2 calculate the delta volley and response (given the 4 values in #1)
    #--- peak assignment, sweep assignment, all values (peaks/valleys for pulse/response & deltas)
 #3 do for each peak in a sweep
 #4 do for each sweep in a file
#do pivot long
test <- function(x,y){
  x+y
}
test(11,12)
#SUBSET SWEEP
soi=2
temp_sweep <- data %>%
  select(paste0("sweep_",soi),time)%>%
  data.frame()

#SUBSET ALL SWEEPS
for (i in c(colnames(data))[2:6]) {
  temp_sweep <- data %>%
    select(i,time)%>%
    data.frame()%>%
    head()
  print(temp_sweep)
}

  # getting the threshold 
 baseline<-mean(data$sweep_1)

standard_deviation<-sd(data$sweep_1)

threshold<-baseline + 3.5 * standard_deviation

data$above_threshold <- data$sweep_1 > threshold

# Find the first index where sweep_1 > threshold
first_above_threshold_index <- which(data$sweep_1 > threshold)[1]
first_above_threshold <- data$time[first_above_threshold_index]


#Finding the end of the Pulse, increasing the numbers will make it 
# more sensitive to finding the end of the pulse, pulse will be shorter

upper_base_limit<-baseline+0.8*standard_deviation

lower_base_limit<-baseline-0.8*standard_deviation

data$baseline_range <- (data$sweep_1 > lower_base_limit) & (data$sweep_1 < upper_base_limit)


## Find the starting index where above_threshold is TRUE
start_index <- which(data$above_threshold == TRUE)[1]

# Initialize a variable to store end_pulse
end_pulse <- NA

# Loop to find 7 consecutive TRUE values in baseline_range
for (i in start_index:(nrow(data) - 6)) {
  if (all(data$baseline_range[i:(i + 6)])) {
    end_pulse <- data$time[i + 6]
    break
  }
}

# Print the end_pulse value
print(paste("End pulse:", end_pulse))

data$derivative_smoothed_sweep1 <- numerical_derivative(data$time, data$smoothed_sweep1)


# Create a new column for when the derivative is close to 0.
data$deriv_0 <- (data$derivative_smoothed_sweep1 > -1) & (data$derivative_smoothed_sweep1 < 1)



# creates a new column dictating which times are part of the first pulse and which are not
data$part_of_pulse <- (data$time > first_above_threshold) & (data$time < end_pulse)

pulse_only_data<-data %>%
  filter(part_of_pulse==TRUE)

data_threshold<-data %>%
  filter(above_threshold==TRUE) 
  
# new graph with first pulse colored blue
ggplot(data,aes(x=time, y=sweep_1,color=part_of_pulse)) +
  geom_line() +
  xlim(0.05,0.2) +
  ylim(-0.04,0.2)
  
## Getting Fiber Volley and Response

# Smooth the data using loess with adjusted span and add to new column
loess_fit <- loess(sweep_1 ~ time, data = data, span = 0.001)  # Adjust span value here
data <- data %>%
  mutate(smoothed_sweep1 = fitted(loess_fit))

ggplot(data,aes(x=time, y=smoothed_sweep1,color=part_of_pulse)) +
  geom_point() +
  xlim(0.058,0.071) +
  ylim(-0.04,0.05)


# Function to calculate numerical derivatives
numerical_derivative <- function(x, y) {
  n <- length(y)
  
  # Check if the lengths of x and y are the same
  if (length(x) != n) {
    stop("Length of x and y must be the same.")
  }
  
  # Initialize the derivative vector
  dy_dx <- numeric(n)
  
  # Forward difference for the first point
  dy_dx[1] <- (y[2] - y[1]) / (x[2] - x[1])
  
  # Central differences for internal points
  for (i in 2:(n - 1)) {
    dy_dx[i] <- (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
  }
  
  # Backward difference for the last point
  dy_dx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  
  return(dy_dx)
}

# Calculate the derivatives of smoothed_pulse1 with respect to time
data$derivative_smoothed_sweep1 <- numerical_derivative(data$time, data$smoothed_sweep1)


# Create a new column for when the derivative is close to 0.
data$deriv_0 <- (data$derivative_smoothed_sweep1 > -1) & (data$derivative_smoothed_sweep1 < 1)



#do the same thing for pulse_only_data
pulse_only_data$derivative_smoothed_sweep1 <- numerical_derivative(data$time, data$smoothed_sweep1)


pulse_only_data$deriv_0 <- (data$derivative_smoothed_sweep1 > -1) & (data$derivative_smoothed_sweep1 < 1)



ggplot(data,aes(x=time, y=smoothed_sweep1,color=deriv_0)) +
  geom_point() +
  xlim(0.058,0.071) +
  ylim(-0.04,0.05)





#Put deriv_0=TRUE point on unsmoothed graphics
ggplot(data,aes(x=time, y=sweep_1,color=deriv_0)) +
  geom_point() +
  xlim(0.058,0.071) +
  ylim(-0.04,0.05)


#For every event (span of trues) find the middle (look inside the pulse only)

fiber_volley_and_pulse<-pulse_only_data%>%
  filter(deriv_0==TRUE)

#If two time entries in the table are close to each other, find the average sweep_1 value for those entries

#Split up table into sections of ten entries each, if the section has deriv_0==TRUE values, average their sweep_1 values

## average all deriv_0==TRUE entries in first 10 entries in pulse_only_data

# Create groups of 10 entries
data_grouped <- pulse_only_data %>%
  mutate(group = ceiling(row_number() / 10))  # Create a grouping variable

# Within each group, filter for deriv_0 == TRUE and find the value furthest from baseline
final_data <- data_grouped %>%
  group_by(group) %>%
  filter(deriv_0 == TRUE) %>%
  mutate(abs_diff_from_baseline = abs(sweep_1 - baseline)) %>%
  filter(abs_diff_from_baseline == max(abs_diff_from_baseline)) %>%
  ungroup() %>%
  select(-group, -abs_diff_from_baseline)  # Remove the grouping and temporary columns

# Print the resulting data frame
print(pulse_only_data)

final_data_final <- head(final_data, 4)

print(final_data_final)

##########################################





ggplot(data, aes(x = time, y = sweep_1, color = above_threshold)) +
  geom_line(data = subset(data, above_threshold == TRUE), color = "blue") +  # Line when above_threshold is TRUE
  geom_line(data = subset(data, above_threshold == FALSE), color = "gray") +  # Line when above_threshold is FALSE
  xlim(0.05, 0.2) +
  ylim(-0.04, 1) +
  scale_color_manual(values = c("blue", "gray"), labels = c("TRUE", "FALSE")) +  # Color scale
  labs(color = "above_threshold") +  # Legend label
  theme_minimal()

ggplot(data,aes(x=time, y=sweep_1, color=above_threshold))+
  geom_point()

#create pluse length variable

#after certain number of points within certain number of standard deviations of baseline, pulse is over

#find volley and response peaks and valleys by smoothing data, then find when deriv=0




ggplot(data, aes(x = time, y = sweep_1)) +
  geom_line(aes(color = above_threshold)) +
  scale_color_manual(values = c("blue", "gray"), labels = c("TRUE", "FALSE")) +
  labs(color = "above_threshold") +
  theme_minimal() +
  xlim(0.05, 0.2) +
  ylim(-0.04, 1) +
  geom_vline(xintercept = data$time[transition_indices] + 0.01, color = "blue", linetype = "dashed")

