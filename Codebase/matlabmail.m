function recipient = matlabmail(recipient, subject, message, attachments) %sender, psswd
% MATLABMAIL Send an email from a predefined gmail account.
%
% MATLABMAIL( recipient, message, subject )
%ab
%   sends the character string stored in 'message' with subjectline 'subject'
%   to the address in 'recipient'.
%   This requires that the sending address is a GMAIL email account.
%
% MATLABMAIL( recipient, message, subject, sender, passwd ) 
%
%   avoids using the stored credentials.
%
% Note: Authentication failed when my gmail account had 2-step verification enabled.
%
% Example:
%
%  There's no example because we don't know your email address!
%  Try to adapt the following:
%
%    pause(60*(1+randi(5))); matlabmail('root@localhost.com', 'done pausing', 'command complete');
%
% See also SENDMAIL


% if nargin<4
%     sender = 'zgimastl@gmx.com';
%     psswd = 'Keko(35)TkO';
% end

sender = 'zgimastl@gmx.com';
psswd = 'Keko(35)TkO';

setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','mail.gmx.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(recipient, subject, message, attachments);


