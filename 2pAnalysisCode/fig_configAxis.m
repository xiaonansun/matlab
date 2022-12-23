function  ax = fig_configAxis(old_axis)

ax = old_axis;

ax.Box = 'off';
ax.Color = 'none';
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
cnt = 1;
for i=1:length(ax.Children)
   if strcmpi(ax.Children(i).Type,'line')
       ax.Children(i).LineWidth = 1;
       if cnt == 1
           ax.Children(i).Color = 'k';
           cnt = cnt + 1;
       end
   end 
end


% Separate x and y axis
% clf
% ax=axes('Position',[.12 .12 .85 .85]);
% hold on
% t=linspace(0,4*pi,100);
% plot(ax,sin(t))
% ax_x=ax;
% ax_x.Position=ax.Position;
% ax_x.Position(2)=.07;
% ax_x.Position(4)=.00001;
% ax_y=ax;
% ax_y.Position=ax.Position;
% ax_y.Position(1)=.07;
% ax_y.Position(3)=.00001;
% ax_x.FontSize = ax.FontSize;
% ax_y.FontSize = ax.FontSize;
% linkprop([ax ax_x],'FontSize');
% linkaxes([ax ax_x],'x')
% linkprop([ax ax_y],'FontSize');
% linkaxes([ax ax_y],'y')
% ax.XColor='none';
% ax.YColor='none';
% ax.XTick=[];
% ax.YTick=[];

