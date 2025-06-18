# VandaMicrobes

Spatial statistical analyses of microbial pinnacles from Antarctic Lakes.

## Project Workflow 

1. [Video data] of benthic microbial mats from Lake Vanda, Antarctica was taken during a 2013 field season.
2. The video data was used to create photogrammetric models using Agisoft Metashape from 19 of the video sample sites. [Model data is contained here.]
3. Pinnacles were extracted from the microbial mat models using [LidarViewer](https://pubs.geoscienceworld.org/gsa/geosphere/article/9/3/546/132601/Point-based-computing-on-scanned-terrain-with) (automatic seperation methods failed due to noise in the dataset).
4. The height and volume of individual pinnacles was extracted and examined. Results and scripts for individual pinnacle analyses are contained within the BLANK folder. 
5. The spatial distribution of the pinnacles was examined by first converting the pinnacles into a point pattern using the maximum height of each pinnacle, then performing spatial analyses on the point pattern. Results for spatial patterning analyses are contained within the BLANK folder.
