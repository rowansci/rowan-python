
from pydantic import BaseModel

from .utils import api_client


class Organization(BaseModel):
    name: str
    weekly_credits: float | None = None
    credits: float | None = None


class OrganizationRole(BaseModel):
    name: str


class SubscriptionPlan(BaseModel):
    name: str


class IndividualSubscription(BaseModel):
    subscription_plan: SubscriptionPlan


class User(BaseModel):
    uuid: str
    username: str
    email: str
    firstname: str | None = None
    lastname: str | None = None
    weekly_credits: float | None = None
    credits: float | None = None
    billing_name: str | None = None
    billing_address: str | None = None
    credit_balance_warning: float | None = None
    organization: Organization | None = None
    organization_role: OrganizationRole | None = None
    individual_subscription: IndividualSubscription | None = None

    @property
    def name(self) -> str:
        return f"{self.firstname or ''} {self.lastname or ''}".strip()

    @property
    def organization_name(self) -> str | None:
        return self.organization.name if self.organization else None

    @property
    def organization_weekly_credits(self) -> float | None:
        return self.organization.weekly_credits if self.organization else None

    @property
    def organization_credits(self) -> float | None:
        return self.organization.credits if self.organization else None

    @property
    def organization_role_name(self) -> str | None:
        return self.organization_role.name if self.organization_role else None

    @property
    def subscription_plan_name(self) -> str | None:
        return (
            self.individual_subscription.subscription_plan.name
            if self.individual_subscription
            else None
        )

    def credits_available_string(self) -> str:
        individual_credits = f"Weekly Credits: {self.weekly_credits}\nCredits: {self.credits}"
        if self.organization is not None:
            return (
                individual_credits
                + f"\nOrganization Weekly Credits: {self.organization_weekly_credits}\n"
                f"Organization Credits: {self.organization_credits}"
            )

        return individual_credits


def whoami() -> User:
    """
    Returns the current user
    """

    with api_client() as client:
        response = client.get("/user/me")
        response.raise_for_status()

    return User(**response.json())
