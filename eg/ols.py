from apop import *
apop_text_to_db("data", "d", 0, 1, None)
data = apop_query_to_data("select * from d")
est  = apop_estimate(data, avar.apop_ols)
print est

#And one last example of using a method of the apop_data object:

print est.parameters.transpose()
