C     COMMON BLOCKS FOR LAMELL MC PROGRAM                         
C                                                                  
      PARAMETER (mion=600,mxspec=10)                                
      CHARACTER*16 title                                              
      REAL ran2
      COMMON/ARR1/ esa(mion,mion),ei(mion)
     *, caver(mxspec),hmik(mxspec,mxspec)
     *, cwi(25,mxspec,mxspec),dum(25)
     *, x6(mion),y6(mion),z6(mion),chv(mion),dp(mion)
     *, rw2(mion),rwi(mion),hc2v(mion,mxspec),hion(mxspec,5),ispc(mion)
     *, pi,avno,bk,eps,epsx,ecf,ech,abeta,dtemp
     *, tx6,ty6,tz6,xww1,du
     *, pwiwal,pwimid,pwiwalv,pwimidv,pcollav,pcollv
     *, box,box2,box2i
     *, ny1,ny2,ny3,my3,nspec,ispec,npart,nwins,nwint  
     *, nfix,il,isos,islu,ink
     *, lll,mmm,jjj,iii,kkk
