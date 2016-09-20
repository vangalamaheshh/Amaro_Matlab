function new_manual_review_station

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Input parameter should be a structure.'); end

P=impose_default_value(P,'reviewing_workstation_IP','nobelium.broadinstitute.org');
P=impose_default_value(P,'IGV_listening_port',60151);
P=impose_default_value(P,'IGV_window_width',2000);


import java.util.*;
import java.io.*;
import java.net.Socket;


[a b] = system('/broad/tools/Linux/x86_64/pkgs/igv/igv &');
if a~=0, error('Problem starting IGV'); end

try
  socket = Socket(P.reviewing_workstation_IP,P.IGV_listening_port);
  IGVout = PrintWriter(socket.getOutputStream(),true);
  IGVin = BufferedReader(InputStreamReader(socket.getInputStream()));
catch me
  fprintf('Failed to establish communication with IGV.  Error was:\n');
  disp(me.message);
  fprintf('Cannot use automatic IGV control\n');
  P.use_IGV_control = false;
end

chrstring = 'chr17';
left = 10e6;
right = left+2000;
IGVout.println(['goto ' chrstring ':' num2str(round(left)) '-' num2str(round(right))]);


