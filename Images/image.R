# load GenVisR and create quality check images
library(GenVisR)

# waterfall plot
png(filename='./Images/waterfall.png', height=10, width=16, units='in', res=150)
waterfall(brcaMAF, main.recurrence_cutoff=.1)
dev.off()