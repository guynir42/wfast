import util.text.sa;

str = ['mex ' sa(getenv('MAT'), '+util/+img/mexCutout.cpp') ' -outdir ' sa(getenv('MAT'), '+util/+img')];

eval(str);