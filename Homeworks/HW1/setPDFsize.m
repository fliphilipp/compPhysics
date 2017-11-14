f = gcf;
f.PaperPositionMode = 'auto';
f_pos = f.PaperPosition;
f.PaperSize = [f_pos(3) f_pos(4)];
