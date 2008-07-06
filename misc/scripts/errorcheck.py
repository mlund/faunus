import sys

tol=1
refval=0.5
val=0.51
exe="Widom"
test="Excess"

err=abs((refval-val)/refval)*100
if (err<=tol): rc="Passed"
else: rc="Failed"

print "%-12s %-12s %-6.1f %-8s" % (exe, test, err, rc)

