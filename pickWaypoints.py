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
# make a class for picking lines:
class LineBuilder:
    def __init__(self,ax,map):
        self.ax = ax
        self.map = map
        self.line,=self.ax.plot(0,0,'r',lw=1.5,zorder=19)
        self.xs = np.array([])
        self.ys = np.array([])
        self.cid = self.ax.figure.canvas.mpl_connect('button_press_event', self)
        self.cids = self.ax.figure.canvas.mpl_connect('pick_event', self.stationpress)
        self.kpid = self.ax.figure.canvas.mpl_connect('key_press_event', self.keypress)
        self.markers,=self.ax.plot(0,0,'rd',picker=5,zorder=20)
        self.edit = False
        self.delete = False
        self.picksta = False
        self.editind=None
        self.name=[]
        self.picked=False
        self.txt=[self.ax.text(0,0,'')]
        # load the stations in case we want to use those...
        self.stations=pickle.load(open('Stations.pickle','r'))    
        self.staline, = map.plot(self.stations['lon'],self.stations['lat'],'d',
                                 color=[0.4,1.,0.4],latlon=True,picker=5,zorder=2)
                

    def keypress(self,event):
        if event.inaxes!=self.ax: return
        print 'keypress', event.key
        if (event.key=='cmd+z') | (event.key=='super+z'):
            print 'Undo!'
            if len(self.xs)>0:
                self.xs=self.xs[:-1]
                self.ys=self.ys[:-1]
                self.plot()
        if (event.key=='e'):
            if not (self.edit):
                print 'Edit!  Press e again to stop'
                self.ax.figure.canvas.mpl_disconnect(self.cid)
                self.ax.figure.canvas.mpl_disconnect(self.cids)
                self.cid = self.ax.figure.canvas.mpl_connect('pick_event', self.editpress)
                self.cid2 = self.ax.figure.canvas.mpl_connect('motion_notify_event', self.editmove)
                self.cid3 = self.ax.figure.canvas.mpl_connect('button_release_event', self.editrelease)
                self.edit=True
            else:
                self.edit=False
                self.ax.figure.canvas.mpl_disconnect(self.cid)
                self.ax.figure.canvas.mpl_disconnect(self.cid2)
                self.ax.figure.canvas.mpl_disconnect(self.cid3)
                self.cid = self.ax.figure.canvas.mpl_connect('button_press_event', self)
                self.cids = self.ax.figure.canvas.mpl_connect('pick_event', self.stationpress)
                print 'Done edit'
        if (event.key=='d'):
            if not (self.delete):
                print 'Edit!  Press d again to stop'
                self.ax.figure.canvas.mpl_disconnect(self.cid)
                self.ax.figure.canvas.mpl_disconnect(self.cids)
                self.cid = self.ax.figure.canvas.mpl_connect('pick_event', self.deletepress)
                self.delete=True
            else:
                self.delete=False
                self.ax.figure.canvas.mpl_disconnect(self.cid)
                self.cid = self.ax.figure.canvas.mpl_connect('button_press_event', self)
                self.cids = self.ax.figure.canvas.mpl_connect('pick_event', self.stationpress)
                print 'Done edit'
    def editpress(self,event):
        print 'Boo'
        print event.artist
        print event.ind
        self.editind=event.ind
    def deletepress(self,event):
        if event.artist==self.staline:
            return
        print 'Boo'
        print event.artist
        print event.ind
        self.xs=np.delete(self.xs,event.ind)
        self.ys=np.delete(self.ys,event.ind)
        print self.name
        self.name.pop(event.ind)
        print self.name
        print 'len(self.name) %02d'%len(self.name)
        print len(self.xs)
        print len(self.xs)

        self.plot()
        self.report()
    def editmove(self,event):
        if event.inaxes is None: return
        if self.editind is None: return
        x,y=event.xdata, event.ydata
        self.xs[self.editind]=x
        self.ys[self.editind]=y
        self.plot()
    def editrelease(self,event):
        self.editind=None
        self.plot()
        self.report()
    def stationpress(self,event):
        print 'Pick', event
        if event.artist==self.staline:
            print 'Picked a station!'
            self.picked=True
            if len(event.ind)>1:
                ind = event.ind[0]
            else:
                ind=event.ind
            x,y=map(self.stations['lon'][ind],self.stations['lat'][ind])            
            print 'Station: %s'%self.stations['name'][ind]
            self.name.append(self.stations['name'][ind])
            self.xs=np.append(self.xs,x)
            self.ys=np.append(self.ys,y)
            self.plot()
            self.report()
    # default click....
    def __call__(self, event):
        # if this already was "picked" then its a station, and we took care of it in stationpress...
        if self.picked:
            self.picked=False
            return
        print 'click', event
        tb = plt.get_current_fig_manager().toolbar
        print tb.mode
        if event.inaxes!=self.ax: return
        if tb.mode!="": return
        
        self.xs=np.append(self.xs,event.xdata)
        self.ys=np.append(self.ys,event.ydata)
        self.name.append('%02d'%len(self.xs))

        self.plot()
        self.report()
            
        

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
        # remove old
        print self.txt
        for txt in self.txt:
            self.ax.texts.remove(txt)
        self.txt=[]
        #print self.txt
        #print self.name[0]
        #print self.xs[0],self.ys[0]
        
        print 'hi!'
        self.txt=[self.ax.text(self.xs[0],self.ys[0],self.name[0],color='r',fontsize=8,zorder=10000)]
        if len(self.xs)>1:
            self.txt.append([self.ax.text(self.xs[-1],self.ys[-1],self.name[-1],color='r',fontsize=8,zorder=1000)])
            for i in range(1,len(self.xs)-1):
                st = '%02d'%(i+1)
                if len(self.name[i])!=2:
                    tx=self.ax.text(self.xs[i],self.ys[i],' '+self.name[i],color='r',fontsize=8,zorder=1000)
                    self.txt.append(tx)
        #self.line=self.ax.plot(self.xs, self.ys,'m-')
        #self.markers=self.ax.plot(self.xs,self.ys,'d')
        self.ax.figure.canvas.draw()

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
            print('%6s:  %02d %05.2f %c  %03d %05.2f %c  %5.1f nm  %5.1f h %5.1f nm %5.1f h '%(self.name[i],latd[i],latm[i],Eq,lond[i],lonm[i],hem,dist[i],dist[i]/7.,cdist[i],cdist[i]/7.)) 
        print('\n')

def format_coord(x, y):
    return 'x=%.4f, y=%.4f'%(map(x, y, inverse = True))


from plotTopo import plotTopo
fig=plt.figure(figsize=(12,12))
ax = plt.gca()
map=plotTopo(fig)
ax.format_coord = format_coord

#map.pcolormesh()
print 'Ready to pick!'
# pickit!
linebuilder = LineBuilder(ax,map)

plt.show()
