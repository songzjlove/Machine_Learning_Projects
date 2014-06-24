import pygmaps

def make_map(nm1, nm2, data) :
    color = {}
    color['CENTER'] = "#CCCCCC"
    color[nm1] = "#CD5C5C"
    color[nm2] = "#87CEFA"
    print color
    centerpoint = (41.8823, -87.6277)
    mymap = pygmaps.maps(centerpoint[0], centerpoint[1], 12)
    mymap.addpoint(centerpoint[0], centerpoint[1], color['CENTER']); 
    print "added center point (%f, %f)" % (centerpoint[0], centerpoint[1])
    
    print "number of data records: %d" % len(data)
    for r in range(len(data)) :
        datum = data[r]
        if datum[12] in color :
            try :
                lat = float(datum[22])
                lng = float(datum[23])
            except :
                ()
            mymap.addpoint(lat, lng, color[datum[12]])
    mymap.draw('./crimemap.html')
    print len(grid)

