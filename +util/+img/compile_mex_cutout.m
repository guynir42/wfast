import util.text.sa;

str = ['mex ' sa(getenv('WFAST'), '+util/+img/mexCutout.cpp') ' -outdir ' sa(getenv('WFAST'), '+util/+img')];

eval(str);