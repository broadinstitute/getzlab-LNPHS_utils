function x = get_LEGO_colors()
   colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7;];
   t = [repmat([3 3 3 3 2 2 2 2 1 1 1 1],4,1);repmat([6 6 6 6 5 5 5 5 4 4 4 4],4,1)];
   r = colors(:,1); g = colors(:,2); b = colors(:,3);
   x = cat(3,r(t),g(t),b(t));
end
