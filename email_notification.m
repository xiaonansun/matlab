function email_notification(subject_message)
% sends an email notification in the subject line to two of my email addresses

setpref('Internet','E_mail','xisun@cshl.edu');
setpref('Internet','SMTP_Server','email.cshl.edu');
sendmail('xisun@cshl.edu',subject_message);
sendmail('xsun4@northwell.edu',subject_message);
