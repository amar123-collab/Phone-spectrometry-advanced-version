function ex_savefig(f, name)
%EX_SAVEFIG  Save a figure as 300 dpi PNG and as PDF.
try
    exportgraphics(f, [name '.png'], 'Resolution', 300);
    exportgraphics(f, [name '.pdf'], 'ContentType', 'vector');
catch
    print(f, [name '.png'], '-dpng', '-r300');
    print(f, [name '.pdf'], '-dpdf');
end
end
