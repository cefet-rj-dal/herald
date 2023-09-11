
# <img src='https://raw.githubusercontent.com/cefet-rj-dal/herald/master/herald_logo.png' align='centre' height='150' width='129'/> Herald

In streaming time series analysis, it is often possible to observe the occurrence of
a significant change in behavior at a certain point or time interval that characterizes the
occurrence of a change point event. 

Change point detection is related to other prominent
areas of monitoring of concept drifts and time series structural breaks, being the subject
of intensive research and applied in areas such as medical condition monitoring, climate
change detection, and human activity analysis. Online, or real-time, change point detection
algorithms run concurrently with the process they are monitoring, processing each data
point as it becomes available. 

Online change point detection for streaming applications is
a fundamentally different and challenging problem that creates an increasing demand for
advanced machine-learning techniques to deal with concept drifts in large data streams.
Generally, regression-based solutions demand the comparison of fixed parametric models,
and a sufficient number of observations in incoming data, whereas solutions based on
deviations from time series prediction can be affected by nonstationary properties and
slowly changing concepts, which may increase detection lags. 

In particular, this research
addresses the challenges posed by prediction-based approaches, exploring the hypothesis
that online change point detection over the lagged prediction of the inertial component
of a time series contributes to smaller detection lags and robustness to nonstationarity
when compared to current approaches. To confirm it, this research proposes __Herald__.

Herald is a novel prediction approach method for
online change point detection based on the use of time series decomposition for deriving
its characteristic trend component as a proxy of its intrinsic inertia. The derived trend
component is given as input to a time series prediction process with a lagged prediction
horizon. The predicted inertial component empowers early detection of
change points in streaming time series presenting nonstationary properties and concept
drifts.

Experimental results indicate that Herald is able to detect change
points with smaller detection lags and with higher sensitivity when compared to baseline
methods for change point detection. Effects are especially relevant in nonstationary time
series with high volatility.
