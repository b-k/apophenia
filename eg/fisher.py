from apop import *
data = [[15, 13], [18, 19]]
adata = apop_pylist_to_data(data) 

print "display the data"
print adata.matrix

result = apop_test_fisher_exact(adata)
print result

print "Or, just one value:"
print result.get("p value", -1)
