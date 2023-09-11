
# Loading Harbinger -------------------------------------------------------


# install.packages("devtools")
#devtools::install_github("cefet-rj-dal/harbinger")

library(harbinger)

#loading the example database
data(har_examples)

# Harbinger without Nexus -------------------------------------------------




#Using the time series 1
dataset <- har_examples[[1]]
head(dataset)


#ploting serie #1
plot(x = 1:length(dataset$serie), y = dataset$serie)
lines(x = 1:length(dataset$serie), y = dataset$serie)


model <- fbiad() |> # establishing method
  fit(dataset$serie) # fitting the model

# making detections using method
detection <- model |>
  detect(dataset$serie)

# filtering detected events
print(detection |> dplyr::filter(event==TRUE))


# evaluating the detections
evaluation <- evaluate(model, detection$event, dataset$event)
print(evaluation$confMatrix)


# ploting the results
grf <- plot.harbinger(model, dataset$serie, detection, dataset$event)
plot(grf)