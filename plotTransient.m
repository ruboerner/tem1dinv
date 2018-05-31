function h = plotTransient(t, V, err, col, as, ref)
%function h = plotTransient(t, V, err, col, as, ref)

if nargin < 6
    ref = [1 length(t) 1 length(t)];
end

if nargin < 5
    as = false;
end

if nargin < 4
    col = 'b';
end

if nargin < 3
    err = [];
end

Vpos = V;
Vneg = -V;
Vpos(Vpos < 0) = NaN;
Vneg(Vneg < 0) = NaN;

Vpos = log10(Vpos);
Vneg = log10(Vneg);
t = log10(t);

VV = nansum([Vpos Vneg], 2);
maxVV = max(VV);
minVV = min(VV);
meanVV = mean(VV);

if isempty(err)
    if as
        h = plot( ...
            t, Vneg, [col '--'], ...
            t, Vpos, [col '+'], ...
            t(ref(1):ref(2)), -5/2 * (t(ref(1):ref(2)) - t(ref(3))) + VV(ref(3)), 'r--', ...
            t(ref(2):end), -1/2   * (t(ref(2):end) - t(ref(4))) + VV(ref(4)), 'b--');
    else
        h = plot( ...
            t, Vneg, [col '--'], ...
            t, Vpos, [col '+']);
    end
else
    err = 1/log(10) * err ./ abs(V);
    hold on;
    h1 = errorbar(t, Vneg, err, col);
    set(h1, 'Marker', 'o', 'LineStyle', 'none');
    h2 = errorbar(t, Vpos, err, col);
    set(h2, 'Marker', '+', 'LineStyle', 'none');
    if as
        %         plot(t(ref(1):ref(2)), -5/2 * (t(ref(1):ref(2)) - t(ref(3))) + VV(ref(3)), 'r--');
        h3 = plot(t, -5/2 * (t - t(ref(3))) + VV(ref(3)), 'r--');
        h4 = plot(t(ref(2):end), -1/2   * (t(ref(2):end) - t(ref(4))) + VV(ref(4)), 'b--');
        legend([h1(1) h2(1) h3(1) h4(1)], 'neg.', 'pos.', '~t^{-2.5}', '~t^{-0.5}');
    else
        legend([h1(1) h2(1)], 'neg.', 'pos.'); 
    end
    hold off;
    
end

xlim([floor(t(1)) ceil(t(end))]);

xlabel('log(t/s)');
ylabel('log(dB_z/dt / V/(Am^2))');
box('on');
grid('on');

end

function y = nansum(x, dim)
if nargin == 1
    aa = ~isnan(x);
    ff = find(aa == 0);
    if (length(ff) ~= 0) % replace NaNs by zeros
        x(ff) = 0;
    end
    ss = sum(aa);    % find the total number of non-nans in the column
    y = sum(x);      % calculate the sum
    ff = find(ss == 0);
    if (length(ff) ~= 0) % put NaNs where the sum came from a column of NaNs
        y(ff) = NaN;
    end
elseif nargin == 2
    aa = ~isnan(x);
    ff = find(aa == 0);
    if (length(ff) ~= 0) % replace NaNs by zeros
        x(ff) = 0;
    end
    ss = sum(aa, dim);   % find the total number of non-nans in the column
    y = sum(x, dim);     % calculate the sum
    ff = find(ss == 0);
    if (length(ff) ~= 0) % put NaNs where the sum came from a column of NaNs
        y(ff) = NaN;
    end
end
end