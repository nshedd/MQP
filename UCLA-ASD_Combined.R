CTL_forcombined = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_SampleLabels')

new.cluster.ids <- c('CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL',
                     'CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL')
names(new.cluster.ids) <- levels(CTL_forcombined)
CTL_forcombined <- RenameIdents(CTL_forcombined, new.cluster.ids)

saveRDS(CTL_forcombined, '/data/rusers/sheddn/UCLA-ASD/data/CTL_GroupLabel')
