function update_logfile(message,output_log)

disp(message)
txtappend(output_log,[datestr(now,31) ' ' message '\n'])


