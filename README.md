# Schelling_on_GIS
Schelling's Segregation Model run with real world demographic data and maps in NetLogo

Contains the [NetLogo](https://ccl.northwestern.edu/netlogo/) model file `Schelling_on_GIS.nlogo`

The model uses NetLogo's GIS extension. It needs GIS data (e.g., shapefiles) to run. These are provided for the UK Local Authority of Bradford in the corresponding directory `Bradford/`. 

The directory `Data_ZuccLoreRodrPaolSerk2020/` includes exports of full simulation runs done with the model for Bradford exported using NetLogo's [`export-world`](http://ccl.northwestern.edu/netlogo/docs/dict/export-cmds.html). The files are provided analog to the Figures 3, 4, 5, and A1 in the scientific manuskript 

"Exploring the dynamics of neighborhood ethnic segregation with agent-based modelling: an empirical application to Bradford "   
by Carolina V. ZUCCOTTI, Jan LORENZ, Alejandra RODRÍGUEZ SÁNCHEZ, Rocco PAOLILLO, and Selamawit SERKA (2021). 

Any of the csv-files in `Data_ZuccLoreRodrPaolSerk2020/` can be loaded in the Model by clicking the "Import World" button in the model's interface. Then the output can be further investigated via various visualization but also run further, or rerun by shuffling the population before.
