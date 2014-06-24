#!/usr/bin/python
# generate_heatmap(locations) makes a heatmap with Javascript functionality in Google maps
# locations is a list of tuples (latitude, longitude, weight)
# an example file is given in locations.txt


prefix = r"""<!DOCTYPE html>
<html style="height: 100%">
  <head>
<meta charset="utf-8">
<title>heatmap</title>
<script src="https://maps.googleapis.com/maps/api/js?v=3.exp&sensor=false&libraries=visualization"></script>
<style type="text/css">
 html { height: 100% }
 body { height: 100%; margin: 0px; padding: 0px }
</style>  
</head>
<body onload="initialize()" style="height: 100%;">
<script>
    var map, pointarray, heatmap;
    var crimelocations = ["""

suffix = r"""];

    function toggleHeatmap() {
        heatmap.setMap(heatmap.getMap() ? null : map);
    }

    function changeGradient() {
      var gradient = [
        'rgba(0, 255, 255, 0)',
        'rgba(0, 255, 255, 1)',
        'rgba(0, 191, 255, 1)',
        'rgba(0, 127, 255, 1)',
        'rgba(0, 63, 255, 1)',
        'rgba(0, 0, 255, 1)',
        'rgba(0, 0, 223, 1)',
        'rgba(0, 0, 191, 1)',
        'rgba(0, 0, 159, 1)',
        'rgba(0, 0, 127, 1)',
        'rgba(63, 0, 91, 1)',
        'rgba(127, 0, 63, 1)',
        'rgba(191, 0, 31, 1)',
        'rgba(255, 0, 0, 1)'
      ]
      heatmap.setOptions({
        gradient: heatmap.get('gradient') ? null : gradient
      });
    }
    
    function changeRadius() {
      heatmap.setOptions({radius: heatmap.get('radius') ? null : 20});
    }
    
    function changeOpacity() {
      heatmap.setOptions({opacity: heatmap.get('opacity') ? null : 0.2});
    }

    function initialize() {
        var mapOptions = {
            zoom: 12,
            center: new google.maps.LatLng(41.862300, -87.627700),
            mapTypeId: google.maps.MapTypeId.HYBRID
        };

        map = new google.maps.Map(document.getElementById('map'),
            mapOptions);

        pointArray = new google.maps.MVCArray(crimelocations);

        heatmap = new google.maps.visualization.HeatmapLayer({
            data: pointArray
        });

        heatmap.setMap(map);
        heatmap.setOptions({radius: heatmap.get('radius') ? null : 20});
  }
</script>
    <div id="panel">
      <button onclick="toggleHeatmap()">Toggle Heatmap</button>
      <button onclick="changeGradient()">Change gradient</button>
      <button onclick="changeRadius()">Change radius</button>
      <button onclick="changeOpacity()">Change opacity</button>
    </div>
   <div id="map" style="height: 100%; width: 100%;"></div>
</body>
</html>"""


def generate_heatmap(weighted_locs, fn="heatmap.html") :
    with open(fn, 'w') as f :
        f.write("%s\n" % prefix)
        for l in weighted_locs :
            f.write("      {location: new google.maps.LatLng(%s, %s), weight: %s},\n" % (l[0], l[1], l[2]))
        f.write("%s\n" % suffix)


f = open('locations.txt', 'r')
locations = []
for l in f :
    l = l.strip()
    loc = l.split(',')
    locations.append((loc[0], loc[1], loc[2]))


generate_heatmap(locations)
