function fig = tulik_figure(figure_number, show_figures)
%TULIK_FIGURE Select a numbered experiment figure without forcing it visible.
%
% MATLAB's figure(n) makes an existing figure visible. That behavior can
% surface diagnostic windows during long table runs even when the root
% default figure visibility is off. Hidden figures use tags instead of
% numbered selection so they can be reused without being shown.

validateattributes(figure_number, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
validateattributes(show_figures, {'logical'}, {'scalar'});

if show_figures
    fig = figure(figure_number);
    return;
end

figure_tag = sprintf('TULIKHiddenFigure%d', figure_number);
figures = findall(groot, 'Type', 'figure', 'Tag', figure_tag);

if isempty(figures)
    fig = figure('Visible', 'off', 'Tag', figure_tag);
else
    fig = figures(1);
    set(fig, 'Visible', 'off');
    set(groot, 'CurrentFigure', fig);
end
end
