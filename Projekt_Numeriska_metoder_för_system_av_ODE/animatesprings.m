function animatesprings(t, x, y, speed, filename)
%ANIMATESPRINGS  Creates a video visualizing the doublespring system.
%   ANIMATESPRINGS(t,x,y) visualizes the spring dynamics given by the
%   the entries in the vectors t, x and y. The time unit is 1/3 second.
%   The output is written to the file springs.mp4.
%   ANIMATESPRINGS(t,x,y,speed) As above, but sets the time unit to be
%   1/speed. Higher speed = shorter video.
%   ANIMATESPRINGS(t,x,y,speed,filename) As above, but writes output to
%   file <filename>.

if nargin < 5
    filename = 'springs.mp4';
end

if nargin < 4
    speed=3;
end

warnState = warning;
warning('off', 'all');

N = length(t);
T = t(end);

% Spring coordinates

m=16;
xw=T/20;
xx=zeros(m+2,1);
xx(3:2:end-2)=-xw;
xx(4:2:end-1)=xw;
yy = linspace(0,1,m+2)';

% Convert differential x,y to actual x,y

xm = min(x); ym = max(y);
L = ym-xm+1;
x = x-L;
y = y-2*L;

ymin = min(y);

fig = figure('Visible','off');   % Don't show plots, faster

% Setup video writer

v = VideoWriter(filename, 'MPEG-4');
fps = 24;
step = round(speed*length(t)/(t(end)*fps)); % Make 1 input sec be 1/speed output secs
v.FrameRate = fps;

% Create video

open(v);

totfr=length(1:step:N);
cnt=1;

fprintf("\n");

for i = 1:step:N
    msg = sprintf("Skapar bildruta %d av %d", cnt, totfr);
    fprintf("%s", msg);

    cnt=cnt+1;
    tt = t(1:i)-t(i)+T;

    plot(T, x(i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); hold on
    plot(T, y(i), 'co', 'MarkerSize', 16, 'MarkerFaceColor', 'c');
    plot(T+xx, x(i)*yy,'-k');
    plot(T+xx, y(i)+(x(i)-y(i))*yy,'-k');
    plot(tt,x(1:i),tt,y(1:i));
    hold off

    axis([-1 T+xw*2 ymin-1 0]);
    axis off
    frame = getframe(fig);
    writeVideo(v, frame);
    fprintf(repmat('\b', 1, strlength(msg)));
end
close(v);

warning(warnState);

fprintf("Skapat videofilen %s\n", filename);

end
