function varargout=init_grad(gradname,type)

grad.head.name          =gradname;
grad.head.points        =0;
grad.head.resolution    =4e-6;  
grad.head.strength      ='';
grad.head.skip          =0;
grad.head.duration      =0;



grad.name               =gradname;
grad.amp                =0;
grad.duration           =4e-6;
grad.tramp              =0;
grad.plateauduration    =0;
grad.resolution         =4e-6;
grad.skip               =0;
grad.m0                 =0;

grad.path               ='/home/Zhiyongi/vnmrsys/shapelib';
grad.gamaHz             =4257.4;
grad.gmax               =100;                                           % [Gauss/cm]
grad.trise              =0.00015;                                       % [s]
grad.slewrate=0.8*grad.gmax/grad.trise;                        % [Gauss/cm/s]


varargout{1}=grad;

