function number = saveReportImage(handler,dir,lab,number,format,fsav)

    if fsav
        saveas(handler,sprintf('%s/%s_%s',dir,lab,num2str(number)),format) 
        number = number + 1;
    end
end