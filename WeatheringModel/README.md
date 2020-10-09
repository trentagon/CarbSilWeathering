## Using the model

The model is fully contained in the weatheringmodel.py file. There are two
example scripts that show how to use the model, simple_example.py and
through_time_example.py. 

The model is configured through the ModelInputs class by setting the various
attributes. The description for each attribute are given in the papers
describing the model. 

To run the model for the modern Earth, simply import the weatheringmodel.py
file, create an input object, then and call the weathering model. For example:

```python
import weatheringmodel
inputs = weatheringmodel.ModelInputs()
results = weatheringmodel.runWeatheringModel(inputs)
```

For details on the inputs and outputs to the model, see the
runWeatheringModel() function in weatheringmodel.py.


