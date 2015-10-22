import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
from mpl_toolkits.basemap import Basemap
import numpy as np
#import seawater as sw
# load the map
#map = Basemap(projection='merc',llcrnrlon=-130.5,llcrnrlat=52.5,urcrnrlon=-127.5,urcrnrlat=54.2,resolution='f')
#pickle.dump(map,open('map.pickle','wb'),-1) 
print 'Loading map'
map = pickle.load(open('map.pickle','rb'))
# load the bathymetry
print 'Loading bathy'
fine=pickle.load(open('coarsebathy.pickle','rb'))
clims = [0,800]

fig=plt.figure(figsize=(12,12))
# make a class for picking lines:
class LineBuilder:
    def __init__(self,ax,map):
        self.ax = ax
        self.map = map
        self.line,=self.ax.plot(0,0,'r')
        self.xs = np.array([])
        self.ys = np.array([])
        self.cid = self.ax.figure.canvas.mpl_connect('button_press_event', self)
        self.kpid = self.ax.figure.canvas.mpl_connect('key_press_event', self.keypress)
        self.markers,=self.ax.plot(0,0,'rd')

    def keypress(self,event):
        if event.inaxes!=self.ax: return
        print 'keypress', event.key
        if (event.key=='cmd+z') | (event.key=='super+z'):
            print 'Undo!'
            if len(self.xs)>0:
                self.xs=self.xs[:-1]
                self.ys=self.ys[:-1]
                self.plot()

    def plot(self):
        # remove old
        print self.markers
#        self.ax.lines.remove(self.line)
#        self.ax.lines.remove(self.markers)
        #plot new
        self.line.set_xdata(self.xs); 
        self.line.set_ydata(self.ys); 
        self.markers.set_xdata(self.xs); 
        self.markers.set_ydata(self.ys); 
        #self.line=self.ax.plot(self.xs, self.ys,'m-')
        #self.markers=self.ax.plot(self.xs,self.ys,'d')
        self.ax.figure.canvas.draw()
        self.report()

    def report(self):
        # make a chart of lat, lon, dist, time...
        lon,lat = map(self.xs,self.ys,inverse=True)
        print lon,lat
        if np.median(lon)<0:
            lon = - lon
            hem = 'W'
        else:
            hem = 'E'
        if np.median(lat)<0:
            lat = -lat
            Eq  = 'S'
        else:
            Eq = 'N'
        lond = np.floor(lon)
        lonm = (lon -lond)*60.
        latd = np.floor(lat)
        latm = (lat-latd)*60
        #dist = sw.dist(lat,lon,units='nm')
        dist = (np.diff(lat)*60)**2 + (np.diff(lon)*60*np.cos(53.5*np.pi/180.))**2
        dist = np.sqrt(dist)
        dist = np.append(0.,dist)
        cdist=np.cumsum(dist)
        for i in range(len(lon)):
            print('%02d:  %02d %05.2f %c  %03d %05.2f %c  %5.1f nm  %5.1f h %5.1f nm %5.1f h '%(i+1,latd[i],latm[i],Eq,lond[i],lonm[i],hem,dist[i],dist[i]/7.,cdist[i],cdist[i]/7.)) 
        print('\n')
    def __call__(self, event):

        print 'click', event
        tb = plt.get_current_fig_manager().toolbar
        print tb.mode
        if event.inaxes!=self.ax: return
        if tb.mode!="": return
        
        self.xs=np.append(self.xs,event.xdata)
        self.ys=np.append(self.ys,event.ydata)
        self.plot()

ax = plt.gca()
def format_coord(x, y):
    return 'x=%.4f, y=%.4f'%(map(x, y, inverse = True))
ax.format_coord = format_coord



print fine.keys()
map.drawcoastlines()
X,Y = np.meshgrid(fine['lon'][:-1],fine['lat'][:-1])
im=map.pcolormesh(X,Y,fine['dep'],rasterized=True,latlon=True,cmap='ocean_r')
plt.clim(clims)
map.contour(X,Y,fine['dep'],levels=[50,100,200,300],colors='0.')
parallels = np.arange(0,100,30./60.)
# labels = [left,right,top,bottom]
map.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,30./60.)
map.drawmeridians(meridians,labels=[True,False,False,True])
map.fillcontinents()
#x,y = map(-130.5,52.5)
#x2,y2 = map(-129,53.5)
#ax.set_xlim(x,x2)
#ax.set_ylim(y,y2)
map.drawparallels(parallels,labels=[False,True,True,False])
map.drawmeridians(meridians,labels=[True,False,False,True])
#map.pcolormesh()
print 'Ready to pick!'
# pickit!
linebuilder = LineBuilder(ax,map)

plt.colorbar(im,shrink=0.5)
plt.show()
