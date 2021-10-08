"""
This file is intended to be used to email job competions/starts and other things. 

Recommened to be used through a bash script that calls this method via python3

"""


class Caclualations_and_jobs():
    def __init__(self) -> None:
        from datetime import datetime
        import os
        self.hname = os.uname()
        print(self.hname)


        # datetime object containing current date and time
        time_now = datetime.now()
        self.dt_finished_string = time_now.strftime("%d/%m/%Y %H:%M:%S")  # date adn time that the script was called (The run finished).

    def message_text(self, ):
        text = f"Hello ,\n\n\
    The calcualtion has completed at {self.dt_now_string} on host machine {self.hname}.\n\nThank you\n\nThis is an automated message\n\n\
    "

    def email(self, attachments = None):
        """
        send_to_actual: If True it will send to all the actual recipients.
        """
        import email, smtplib, ssl
        from email import encoders
        from email.mime.base import MIMEBase
        from email.mime.multipart import MIMEMultipart
        from email.mime.text import MIMEText

        port = 587  # For starttls
        smtp_server = "smtp.live.com"
        sender_email = statusreport_eranjan@outlook.com
        password = my_password
        context = ssl.create_default_context()
        with smtplib.SMTP(smtp_server, port) as server:
            server.ehlo()  # Can be omitted
            server.starttls(context=context)
            server.ehlo()  # Can be omitted
            server.login(sender_email, password)

            # Some inits
            email_names = []
            email_cc = []
            email_emails = []
            email_salutation = []
            email_institution = []

            # Getting all data
            if send_to_actual:
                for i in range(num_rows):
                    email_emails.append(df.iloc[i]["to"])
                    email_cc.append(df.iloc[i]["cc"])
                    email_names.append(df.iloc[i]["name"])
                    email_salutation.append(df.iloc[i]["salutation"])
                    email_institution.append(df.iloc[i]["institution"])
            else:
                for i in range(num_rows):
                    email_emails.append("maduhaseeban@outlook.com")
                    email_cc.append("ebk_era@yahoo.com")
                    email_names.append(df.iloc[i]["name"].strip())
                    email_salutation.append(df.iloc[i]["salutation"])
                    email_institution.append(df.iloc[i]["institution"].strip())

            # Here we start the emailing process
            for n in range(0,len(email_names)):
                rcpts = []
                name = email_names[n]
                receiver_email = email_emails[n]
                salutation = email_salutation[n]
                rcpts.append(receiver_email)
                # print(f"Sending email to {name}  - ({n+1}/{len(email_names)})")
                body  = message_text(salutation, name, email_institution[n])
                if type(email_cc[n]) != float:
                    print(f"Emailing: {name};\t@:{receiver_email}\tCC @: {email_cc[n]}\t- ({n+1}/{len(email_names)})\t Actual email: {send_to_actual}")
                else:
                    print(f"Emailing: {name};\t@:{receiver_email}\t\t- ({n+1}/{len(email_names)})\t Actual email: {send_to_actual}")

                # Create a multipart message and set headers
                message = MIMEMultipart()
                message["From"] = "Status Report"
                message["To"] = receiver_email
                if type(email_cc[n]) != float:
                    message["CC"] = email_cc[n].strip()  # Recommended for mass emails
                    rcpts.extend(email_cc[n].strip().split(","))
                # message["Bcc"] = receiver_email  # Recommended for mass emails
                message["Subject"] = f"Applying for internal medicine residency position at {email_institution[n]}"

                # Open PDF file in binary mode
                # print(attachments)
                if attachments:
                    try:
                        for filename in attachments:
                            # karapu wenasa thamai attachment kiyana eka attachments karala menna me for loop eka aluthn ekathu kara
                            with open(filename, "rb") as attachment:
                                filename = filename.strip("../")
                                print(f"Attaching: {filename}")
                                # Add file as application/octet-stream
                                # Email client can usually download this automatically as attachment
                                part = MIMEBase("application", "octet-stream")
                                part.set_payload(attachment.read())
                            # Encode file in ASCII characters to send by email    
                            encoders.encode_base64(part)
                            # Add header as key/value pair to attachment part
                            part.add_header("Content-Disposition",f"attachment", filename=filename)
                            # Add attachment to message and convert message to string
                            message.attach(part)
                    except:
                        print(f"Could not find attachment: {filename}, ignoring recipient: {name} and stopping sending emails")
                        break

                # Here goes all the email sending part
                # Add body to email
                message.attach(MIMEText(body, _charset="UTF-8"))
                # print(message)
                text = message.as_string()
                # print(rcpts)  # For debugging
                # print(text)  # For debugging
                # print("E mail sent!")  # Leaving this here for future debugging purposes
                server.sendmail(sender_email, rcpts, text)

