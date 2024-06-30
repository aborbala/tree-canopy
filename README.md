# Individual Tree Crowns Segmentation Using Deep Learning Methods

## Data Sources

### DOP 2020
- **Description**: Digital Orthophoto (DOP) provides high-resolution (20cm) aerial imagery that serves as a base layer for various geospatial analyses. Tiles: 2km*2km.  Bildflug 01., 08., 12. und 16. August 2020 
- **Source**: [[DOP Source Link](https://daten.berlin.de/datensaetze/digitale-farbige-orthophotos-2020-dop20rgb-wms)]

### LAS
- **Description**: LAS files contain LiDAR point cloud data, which includes detailed 3D information about the Earth's surface, vegetation, and structures. Airborne Laserscanning Flug 24.02.2021, 25.02.2021 und 02.03.2021
- **Source**: [[LAS Source Link](https://fbinter.stadt-berlin.de/fb/berlin/service_intern.jsp?id=a_als@senstadt&type=FEED)]

### Building Footprints (ALKIS)
- **Description**: Building footprints (ALKIS) are vector data representing the outline of buildings. This data is used to eliminate building points from the LiDAR data to avoid confusion with tree crowns.
- **Source**: [[Building Footprints Source Link](https://www.berlin.de/sen/sbw/stadtdaten/geoportal/liegenschaftskataster/)]
- **WFS**: https://fbinter.stadt-berlin.de/fb/wfs/data/senstadt/s_wfs_alkis_gebaeudeflaechen

## 01. Preprocessing Steps

### 1. DOP Slicing
- Slice the DOP into 100m x 100m segments to manage the data efficiently and to match the processing scale of the LiDAR data.

### 2. LAS and Building Footprints Integration
- Combine LAS data with building footprints to eliminate points corresponding to buildings from the LAS data. This step ensures that only vegetation and natural features are analyzed for tree crown extraction.

### 3. LAS Cropping and Tree Crown Extraction
- Using the dimensions of the DOP slices, crop the LAS data to correspond with each 100m x 100m segment.
- Extract tree crowns from each LAS slice based on height and spatial distribution.

### 4. Data Cleaning
- Eliminate polygons with elongated shapes that are unlikely to represent tree crowns. This step helps in refining the dataset to include only tree-like shapes.

## 02. Ranking of Extraction Results
- Manually investigate all slices and rank them based on the quality of tree crown extraction. This ranking helps in selecting the best samples for model training.

## 03. Model Training Using Detectron2
- Utilize Detectron2, a state-of-the-art deep learning library for object detection and segmentation, to train a model on the ranked data. The model learns to identify and segment tree crowns from the input data.

## 04. Evaluation
- Evaluate the trained model on a separate validation set to assess its performance. Metrics such as precision, recall, and F1 score will be used to determine the accuracy and effectiveness of the segmentation model.


